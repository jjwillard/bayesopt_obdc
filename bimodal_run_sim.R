
pacman::p_load(hetGP, ggdist, dplyr, tidyr, purrr, foreach, doParallel, spacefillr, cli)
source("./simulation_functions.R")

#------------------------------------------------------------------------#
# Data generating mechanism
#------------------------------------------------------------------------#
  
# single binary covariate
z_grid <-  expand_grid(z1 = c(0,1))

# yields single value of objective function for specific scenario
objective_function_evaluation <- function(dose, mean_vector_f, covariance_matrix_f,
                                          mvtnorm_scaling_f, base_response_f, constant_value_f){
  
  if(all(is.na(mean_vector_f))){
    constant_value_f
  } else {
    -mvtnorm_scaling_f * mvtnorm::dmvnorm(dose,
                                          mean_vector_f,
                                          covariance_matrix_f) - base_response_f
  }
  
}

constraint_function_evaluation <- function(dose, mean_vector_g, covariance_matrix_g,
                                           mvtnorm_scaling_g, base_response_g, constant_value_g){
  
  if(all(is.na(mean_vector_g))){
    constant_value_g
  } else {
    mvtnorm_scaling_g * mvtnorm::dmvnorm(dose,
                                         mean_vector_g,
                                         covariance_matrix_g) + base_response_g
  }
  
}

get_optimal_dose <- function(dimension_x, x_grid_granularity, z, 
                             mean_vector_f, covariance_matrix_f,
                             mvtnorm_scaling_f, base_response_f, constant_value_f,
                             mean_vector_g, covariance_matrix_g,
                             mvtnorm_scaling_g, base_response_g, constant_value_g,
                             constraint_threshold_g){
  .names_x <- paste0("x", 1:dimension_x)
  .grid_x_single_dimension <- seq(0, 1, x_grid_granularity)
  .grid_x <- map2(.names_x, rep(list(.grid_x_single_dimension), dimension_x),
                  ~ tibble(!!.x := .y)) %>% list_cbind() %>%
    reduce(expand_grid, .name_repair = "minimal") %>%
    set_names(.names_x)
  .grid_x_list <- .grid_x %>% transpose() %>%
    map(~ unlist(.), use.names = TRUE)
  
  #----------------------------------------#
  #----------- CHANGE HERE-----------------#
  #----------------------------------------#
  .values <- .grid_x %>%
    mutate(true_f_at_optimal = map_dbl(.grid_x_list,
                                       objective_function_evaluation,
                                       mean_vector_f = mean_vector_f,
                                       covariance_matrix_f = covariance_matrix_f,
                                       mvtnorm_scaling_f = mvtnorm_scaling_f,
                                       base_response_f = base_response_f,
                                       constant_value_f = constant_value_f),
           true_g_at_optimal = map_dbl(.grid_x_list,
                                       constraint_function_evaluation,
                                       mean_vector_g = mean_vector_g, 
                                       covariance_matrix_g = covariance_matrix_g,
                                       mvtnorm_scaling_g = mvtnorm_scaling_g, 
                                       base_response_g = base_response_g, 
                                       constant_value_g = constant_value_g),
           satisfies_constraint = if_else(true_g_at_optimal <= constraint_threshold_g, 1, 0)) %>%
    filter(satisfies_constraint == 1) %>%
    filter(true_f_at_optimal == min(true_f_at_optimal)) %>%
    select(-satisfies_constraint) %>%
    mutate(constraint_threshold_g = constraint_threshold_g)
  
  .true_optimal_dose <- .values %>% select(num_range("x", 1:20)) %>% transpose() %>% map(~ unlist(.), use.names = TRUE)
  
  .values %>%
    mutate(true_optimal_dose = .true_optimal_dose) %>%
    select(true_optimal_dose, true_f_at_optimal, true_g_at_optimal, constraint_threshold_g) %>%
    bind_cols(z %>% as_tibble() %>% set_names(names(z)), .name_repair = "minimal") #in case more than one optimum
  
}

objective_constraint_parameters <-  z_grid %>%
  #put covariate strata vectors into list form
  mutate(z = z_grid %>%
           transpose() %>%
           map(~ unlist(.), use.names = TRUE),
         
         
         # specify stratum specific distributional information
         # MUST MATCH ORDER OF ROWS OF Z_GRID!!!!
         
         # this is for objective function f
         mean_vector_f = list(c(0.25, 0.75), c(0.75, 0.25)),
         covariance_matrix_f = rep(list(diag(0.1, 2)), nrow(.)),
         mvtnorm_scaling_f = c(1),
         base_response_f = c(0),
         constant_value_f = c(0),
         
         # this is for constraint function g
         mean_vector_g = list(c(0.75, 1.25), c(1.25, 0.75)),
         covariance_matrix_g = rep(list(diag(0.1, 2)), nrow(.)),
         mvtnorm_scaling_g = c(1),
         base_response_g = c(0),
         constant_value_g = c(0),
         constraint_threshold_g = c(0.2,0.2))



true_optimal_dose <- pmap(list(z = objective_constraint_parameters$z,
                               mean_vector_f = objective_constraint_parameters$mean_vector_f,
                               covariance_matrix_f = objective_constraint_parameters$covariance_matrix_f,
                               mvtnorm_scaling_f = objective_constraint_parameters$mvtnorm_scaling_f,
                               base_response_f = objective_constraint_parameters$base_response_f,
                               constant_value_f = objective_constraint_parameters$constant_value_f,
                               mean_vector_g = objective_constraint_parameters$mean_vector_g,
                               covariance_matrix_g = objective_constraint_parameters$covariance_matrix_g,
                               mvtnorm_scaling_g = objective_constraint_parameters$mvtnorm_scaling_g,
                               base_response_g = objective_constraint_parameters$base_response_g,
                               constant_value_g = objective_constraint_parameters$constant_value_g,
                               constraint_threshold_g = objective_constraint_parameters$constraint_threshold_g),
                          get_optimal_dose,
                          dimension_x = 2,
                          x_grid_granularity = 0.25) %>% list_rbind()

distribution_parameters <- objective_constraint_parameters %>%
  left_join(true_optimal_dose %>% select(-constraint_threshold_g), by = names(z_grid))


objective_function <- function(x, z, distribution_parameters){
  .names_x <- paste0("x", 1:length(x))
  .data_x <- as.matrix(map2(.names_x, x, ~ tibble(!!.x := .y)) %>% list_cbind())
  .index_z <- distribution_parameters$z %>%
    map(~ all(length(.) == length(z)) && all(. == z)) %>%
    list_c() %>%
    which()
  
  # Abort if z vector isn't part of data generating mechanism
  if(length(.index_z) == 0L){
    cli_abort(c(message = "Can't find z vector in {.var distribution_parameters$z}. Make sure it's the right length and has valid values.",
                "i" = "z vector: {z}",
                "i" = "x vector: {x}"))
  }
  
  #----------------------------------------#
  #----------- CHANGE HERE-----------------#
  #----------------------------------------#
  objective_function_evaluation(dose = .data_x,
                                mean_vector_f = distribution_parameters$mean_vector_f[[.index_z]],
                                covariance_matrix_f = distribution_parameters$covariance_matrix_f[[.index_z]],
                                mvtnorm_scaling_f = distribution_parameters$mvtnorm_scaling_f[[.index_z]],
                                base_response_f = distribution_parameters$base_response_f[[.index_z]],
                                constant_value_f = distribution_parameters$constant_value_f[[.index_z]])
  
}


constraint_function <- function(x, z, distribution_parameters){
  .names_x <- paste0("x", 1:length(x))
  .data_x <- as.matrix(map2(.names_x, x, ~ tibble(!!.x := .y)) %>% list_cbind())
  .index_z <- distribution_parameters$z %>%
    map(~ all(length(.) == length(z)) && all(. == z)) %>%
    list_c() %>%
    which()
  
  # Abort if z vector isn't part of data generating mechanism
  if(length(.index_z) == 0L){
    cli_abort(c(message = "Can't find z vector in {.var distribution_parameters$z}. Make sure it's the right length and has valid values.",
                "i" = "z vector: {z}",
                "i" = "x vector: {x}"))
  }
  
  #----------------------------------------#
  #----------- CHANGE HERE-----------------#
  #----------------------------------------#
  constraint_function_evaluation(.data_x,
                                 mean_vector_g = distribution_parameters$mean_vector_g[[.index_z]],
                                 covariance_matrix_g = distribution_parameters$covariance_matrix_g[[.index_z]],
                                 mvtnorm_scaling_g = distribution_parameters$mvtnorm_scaling_g[[.index_z]],
                                 base_response_g = distribution_parameters$base_response_g[[.index_z]],
                                 constant_value_g = distribution_parameters$constant_value_g[[.index_z]])
  
}

#------------------------------------------------------------------------#
# Simulation settings
#------------------------------------------------------------------------#

final_filename_rds <- "./Data/bimodal_resp_het_f_g_toxthresh_0.2.RDS"

n_simulations <- 1100

constraint_standard <- tibble(constraint_threshold_g = c(0.2)) # set this separately to illustrate selecting only one for standard case
constraint_personalized <- distribution_parameters %>% 
  select(num_range("z", 1:20), constraint_threshold_g) 

design_parameters <- tibble(initial_n_doses = c(1, 2),
                            reps_each_x =  c(4, 2),
                            algorithm = c("standard", "personalized"),
                            constraint_threshold_g = list(constraint_standard, constraint_personalized),
                            z_enrollment_type = c("random", "balanced"),
                            maximum_sample_size = 80,
                            x_generation_type = c("dose_escalation", "dose_escalation"),
                            x_grid_granularity = c(0.25),
                            acquisition_function = "CEI",
                            constraint_probability_threshold = c(0.9),
                            rho = c(0.25)) 

# looking at 0.79/1/3.77 ses 
additional_parameters <- expand_grid(true_noise_sd_f = c(1.591549),
                                     true_noise_sd_g = c(0.1306423))

design_parameters <- design_parameters %>%
  uncount(nrow(additional_parameters)) %>%
  bind_cols(bind_rows(additional_parameters, additional_parameters)) 


#------------------------------------------------------------------------#
# Standard Algorithm
#------------------------------------------------------------------------#

dose_finding <- safely(dose_finding)

#tictoc::tic()
registerDoParallel(cores = parallel::detectCores())

seed <-  100000:(100000 + n_simulations - 1)

sim_results <- foreach(i = 1:nrow(design_parameters),
                       .errorhandling = "pass",
                       .combine = 'rbind') %do% {
                         
                 foreach(j = 1:n_simulations,
                         .errorhandling = "pass",
                         .combine = 'rbind') %dopar% {
                    
              dose_finding(simulation = j,
                           algorithm = design_parameters$algorithm[i],
                           objective_function = objective_function, 
                           constraint_function = constraint_function,
                           constraint_threshold_g = design_parameters$constraint_threshold_g[[i]],
                           distribution_parameters = distribution_parameters,
                           true_noise_sd_f = design_parameters$true_noise_sd_f[i], 
                           true_noise_sd_g = design_parameters$true_noise_sd_g[i], 
                           initial_n_doses = design_parameters$initial_n_doses[i], 
                           x_grid_granularity = design_parameters$x_grid_granularity[i], 
                           seed = seed[j], 
                           reps_each_x = design_parameters$reps_each_x[i],
                           max_sample_size = design_parameters$maximum_sample_size[i],
                           acquisition_function = design_parameters$acquisition_function[i],
                           x_generation_type = design_parameters$x_generation_type[i],
                           z_enrollment_type = design_parameters$z_enrollment_type[i],
                           af_stop_threshold = -1, #no early stopping
                           constraint_probability_threshold = design_parameters$constraint_probability_threshold[i],
                           rho = design_parameters$rho[i])

                                 }
                       }
                       
stopImplicitCluster()
#tictoc::toc()

saveRDS(sim_results, final_filename_rds)

