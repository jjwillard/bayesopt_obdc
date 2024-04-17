# This looks at a toxicity threshold of 1.5 and 2 for mild/moderate vs severe for the log(aeburden + 0.5) endpoint 

pacman::p_load(hetGP, ggdist, dplyr, tidyr, purrr, foreach, doParallel, spacefillr, cli)
source("./simulation_functions.R")

#------------------------------------------------------------------------#
# Data generating mechanism
#------------------------------------------------------------------------#
  
# single binary covariate
z_grid <-  expand_grid(z1 = c(0,1))

objective_function_evaluation <- function(dose, z){
  d_1 <- dose[1]
  d_2 <- dose[2]

  # Y_f = \beta_0 + \beta_1 d_1 + \beta_2 d_2 + \beta_3 d_1 d_2 + \beta_4 d_1^2 + \beta_5 d_2^2 + \beta_6 d_1^2 d_2^2
  if(z == 0){
    -1.3796 - 4.0794 * d_1 - 0.4827 * d_2 - 4.2250 * d_1 * d_2 + 2.4515 * d_1^2 - 7.5101 * d_2^2 - 1.5605 * d_1^2 * d_2^2

  # severe OSA subgroup  
  } else if (z == 1){
    1.054 - 11.277* d_1 - 8.322 * d_2 - 17.024 * d_1 * d_2  + 8.174 * d_1^2 + 2.338 * d_2^2 + 4.607 * d_1^2 * d_2^2
    
  }
  
}

constraint_function_evaluation <- function(dose){
  d_1 <- dose[1]
  d_2 <- dose[2]

  # Y_g = \beta_0 + \beta_1d_1 + \beta_2d_2 + \beta_3 d_1 d_2 + \beta_4 d_1^2 + \beta_5 d_2^2 + \beta_6 d_1^2 d_2^2
  
  -0.5853044 + 1.8338103 * d_1 + 2.2618823 * d_2 - 4.0500226 * d_1 * d_2 + 
    1.7880889 * d_1^2 + 0.4653626* d_2^2 + 2.9075227 * d_1^2 * d_2^2

}


get_optimal_dose <- function(dimension_x, x_grid_granularity, z,
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
    mutate(true_f_at_optimal = map2_dbl(.grid_x_list,
                                       z,
                                       objective_function_evaluation),
           true_g_at_optimal = map_dbl(.grid_x_list,
                                       constraint_function_evaluation),
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
         constraint_threshold_g = c(1.5,2)) 

true_optimal_dose <- pmap(list(z = objective_constraint_parameters$z,
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
                                z = z)
  
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
  constraint_function_evaluation(.data_x)
  
}

# constraint_function(c(0.25, 0.75), z = c(0), distribution_parameters)
# constraint_function(c(0.5,0.75), z = c(1), distribution_parameters)

#------------------------------------------------------------------------#
# Simulation settings
#------------------------------------------------------------------------#

final_filename_rds <- "rwa_no_early_stopping_toxthresh_1.5_2_rho_0.5.RDS"

n_simulations <- 1100

constraint_standard <- tibble(constraint_threshold_g = c(1.5)) 
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
                            rho = c(0.5),
                            true_noise_sd_f = 7.68,
                            true_noise_sd_g = 1.29) 


#------------------------------------------------------------------------#
# Dose Finding
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

