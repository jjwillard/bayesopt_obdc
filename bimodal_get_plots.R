## Set of functions for plotting both regular/personalized algorithms

pacman::p_load(dplyr, tidyr, purrr, ggplot2, cowplot, foreach, doParallel, latex2exp, cli)

##################### Data manipulation #############################

sim_results_constraint_0.5 <- readRDS("bimodal_resp_het_f_g_toxthresh_0.5.RDS")
results_constraint_0.5 <- sim_results_constraint_0.5[,1] %>%
  list_rbind() 


sim_results_constraint_0.2 <- readRDS("bimodal_resp_het_f_g_toxthresh_0.2.RDS")
results_constraint_0.2 <- sim_results_constraint_0.2[,1] %>%
  list_rbind() 


results <- results_constraint_0.5 %>%
  bind_rows(results_constraint_0.2)

results_no_missing <- results %>% filter(!is.na(f)) %>%
  group_by(algorithm, seed, constraint_threshold_g, true_noise_sd_f, true_noise_sd_g) %>%
  mutate(evaluation = row_number()) %>%
  ungroup()

# need to manually specify iteration number if one strata stopped and other continued
results_no_missing <- results_no_missing %>%
  mutate(iteration2 = case_when(evaluation <=4 ~ 0,
                                evaluation > 4 & evaluation <=8 ~ 1,
                                evaluation > 8 & evaluation <=12 ~ 2,
                                evaluation > 12 & evaluation <=16 ~ 3,
                                evaluation > 16 & evaluation <=20 ~ 4,
                                evaluation > 20 & evaluation <=24 ~ 5,
                                evaluation > 24 & evaluation <=28 ~ 6,
                                evaluation > 28 & evaluation <=32 ~ 7,
                                evaluation > 32 & evaluation <=36 ~ 8,
                                evaluation > 36 & evaluation <=40 ~ 9,
                                evaluation > 40 & evaluation <=44 ~ 10,
                                evaluation > 44 & evaluation <=48 ~ 11,
                                evaluation > 48 & evaluation <=52 ~ 12,
                                evaluation > 52 & evaluation <=56~ 13,
                                evaluation > 56 & evaluation <=60 ~ 14,
                                evaluation > 60 & evaluation <=64 ~ 15,
                                evaluation > 64 & evaluation <=68 ~ 16,
                                evaluation > 68 & evaluation <=72 ~ 17,
                                evaluation > 72 & evaluation <=76 ~ 18,
                                evaluation > 76 & evaluation <=80~ 19))

#--------------- Plotting ---------------------------------------#

tick_labels <- c("1\n8", "5\n24", "10\n44", "15\n64", "19\n80") 

summarized_data <- results_no_missing %>%
  group_by(algorithm, z1, 
           true_noise_sd_f, true_noise_sd_g, iteration2, 
           constraint_threshold_g) %>%
  summarise(mean_rpsel_true_f_at_optimal = mean(rpsel_true_f_at_optimal, na.rm = TRUE),
            mean_rpsel_true_g_at_optimal = mean(rpsel_true_g_at_optimal, na.rm = TRUE),
            mean_rpsel_true_f_at_best_x = mean(rpsel_true_f_at_best_x, na.rm = TRUE),
            mean_rpsel_true_g_at_best_x = mean(rpsel_true_g_at_best_x, na.rm = TRUE),
            mean_f = mean(posterior_mean_f_at_best_x, na.rm = TRUE),
            mean_g = mean(posterior_mean_g_at_best_x, na.rm = TRUE),
            mean_dose_units_from_optimal = mean(dose_units_from_optimal, na.rm = TRUE),
            mean_absolute_deviation_f = mean(abs(posterior_mean_f_at_best_x - true_f_at_best_x), na.rm = TRUE),
            mean_absolute_deviation_g = mean(abs(posterior_mean_g_at_best_x - true_g_at_best_x), na.rm = TRUE),
            .groups = "drop") %>%
         mutate(z1 = factor(z1),
                true_noise_sd_f = factor(true_noise_sd_f),
                true_noise_sd_g = factor(true_noise_sd_g),
                constraint_threshold_g = factor(constraint_threshold_g))
levels(summarized_data$z1) <- c(`0` = TeX("$z_1 = 0$"), `1` = TeX("$z_1 = 1$"))
levels(summarized_data$true_noise_sd_f) <- c(`1.591549` = TeX("$ses_f=1$"))
levels(summarized_data$true_noise_sd_g) <- c(`0.1306423` = TeX("$ses_g=1$"))

levels(summarized_data$constraint_threshold_g) <- c(`0.5` = TeX("$g^{\\dagger}=0.5$"),
                                                    `0.8` = TeX("$g^{\\dagger}=0.8$"))



#----------------------------- Dose units from optimal ---------------------------------#

plot_exp_dose_units <- ggplot(summarized_data %>% filter(true_noise_sd_f == "ses[f] * {\n    phantom() == phantom()\n} * 1", 
                                                        true_noise_sd_g == "ses[g] * {\n    phantom() == phantom()\n} * 1",
                                                        iteration2 >= 1), 
       aes(x = iteration2, y = mean_dose_units_from_optimal, 
           color = as.factor(algorithm),
           linetype = as.factor(algorithm),
           shape = constraint_threshold_g)) + 
  facet_grid(z1 ~ algorithm, labeller = label_parsed, scales = "free") +
  geom_vline(xintercept = 8, linetype = "dashed", alpha = 0.2) +
  geom_point(size = 1.7) +
  geom_line() +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 19),
                     labels = tick_labels) +
  scale_y_continuous(limits = c(NA,NA)) +
  scale_color_manual("algorithm", values = c("#26495c", "#c66b3d")) +
  scale_linetype_manual(name = "algorithm",
                        values = c(11, "dotdash")) +
  scale_shape_manual(parse(text = TeX("$g^{\\dagger}$")),
                     values = c(0,6),
                     labels = c("0.2", "0.5")) +
  theme_bw() +
  guides(color = guide_legend(nrow = 1, order = 1, override.aes=list(linewidth = 1, shape = NA)),
         linetype = guide_legend(nrow = 1, order = 1)) +
  ylab("dose units from optimal") +
  xlab("iteration \n total number of participants") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "top",
        #legend.title = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 7),
        aspect.ratio = 0.75)



#----------------------------------- RPSEL at optimal --------------------------------#
plot_avg_rpsel <- ggplot(summarized_data %>% filter(true_noise_sd_f == "ses[f] * {\n    phantom() == phantom()\n} * 1", 
                                  true_noise_sd_g == "ses[g] * {\n    phantom() == phantom()\n} * 1",
                                  iteration2 >= 1), 
       aes(x = iteration2, y = mean_rpsel_true_f_at_optimal, 
           color = as.factor(algorithm),
           linetype = as.factor(algorithm),
           shape = as.factor(constraint_threshold_g))) + 
  facet_grid(z1 ~ algorithm, labeller = label_parsed, scales = "free") +
  geom_vline(xintercept = 8, linetype = "dashed", alpha = 0.2) +
  geom_point(size = 1.7) + 
  geom_line() +
  scale_x_continuous(breaks = c(1, 5, 10, 15, 19),
                     labels = tick_labels) +
  scale_color_manual("algorithm", values = c("#26495c", "#c66b3d")) +
  scale_linetype_manual(name = "algorithm",
                        values = c(11, "dotdash")) +
  scale_shape_manual(parse(text = TeX("$g^{\\dagger}$")),
                     values = c(0,6),
                     labels = c("0.2", "0.5")) +
  theme_bw() +
  guides(color = guide_legend(nrow = 1, order = 1, override.aes=list(linewidth = 1, shape = NA)),
         linetype = guide_legend(nrow = 1, order = 1)) +
  ylab("rpsel") +
  xlab("iteration \n total number of participants") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "top",
        #legend.title = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 7),
        aspect.ratio = 0.75)



#------------------------ Plot for expected number of toxic doses -----------------------------#
toxic_doses <- results_no_missing %>%
  filter(true_noise_sd_f == 1.591549, true_noise_sd_g == 0.1306423) %>%
  group_by(algorithm, z1, 
           true_noise_sd_f, true_noise_sd_g,
           constraint_threshold_g, seed) %>%
  summarize(sum_evaluated_x_violates_constraint = sum(evaluated_x_violates_constraint, na.rm = TRUE),
            .groups = "keep") %>%
  ungroup(seed) %>%
  summarize(mean_number_parts_toxicity = mean(sum_evaluated_x_violates_constraint),
            .groups = "drop")

toxic_doses_personalized <- toxic_doses %>%
  filter(algorithm == "personalized")

toxic_doses_standard <- toxic_doses %>%
  filter(algorithm == "standard") %>%
  select(-z1) %>%
  group_by(algorithm, true_noise_sd_f, true_noise_sd_g, constraint_threshold_g) %>%
  summarise(mean_number_parts_toxicity = sum(mean_number_parts_toxicity),
            .groups = "drop")

toxic_doses_all <- toxic_doses_personalized %>%
  bind_rows(toxic_doses_standard) %>%
  mutate(x_generation = "dose escalation",
         constraint_threshold_g = factor(constraint_threshold_g))
levels(toxic_doses_all$constraint_threshold_g) <- c(`0.2` = TeX("$g^{\\dagger}=0.2$"),
                                                    `0.5` = TeX("$g^{\\dagger}=0.5$"))

plot_toxic_doses <- ggplot(toxic_doses_all,
                           aes(x = constraint_threshold_g,
                               y = mean_number_parts_toxicity,
                               fill = factor(z1))) +
  facet_grid(. ~ algorithm, labeller = label_parsed) +
  geom_col() +
  scale_fill_manual(values = c("#79AEB2", "#4A6274"),
                    na.value = "#657580",
                    limits = c("0", "1"),
                    labels = unname(TeX(c("$z_1=0$", "$z_1=1$")))) +
  scale_x_discrete(labels = c("0.2", "0.5")) +
  theme_bw() +
  ylab("toxic doses") +
  xlab(parse(text = TeX("$g^{\\dagger}"))) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        #legend.position = "top",
        legend.title = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        aspect.ratio = 1)


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


x_grid <- expand_grid(x1 = seq(0, 1, length.out = 100),
                      x2 = seq(0, 1, length.out = 100))
x_grid <- bind_rows(x_grid, x_grid)

x_list <- x_grid %>%
  purrr::transpose() %>%
  map(tibble::as_tibble_row)

z_grid <- tibble(z1 = c(rep(0, length(x_list)/2),
                        rep(1, length(x_list)/2)))
z_list <- z_grid  %>%
  purrr::transpose() %>%
  map(tibble::as_tibble_row)

f <- map2_dbl(x_list,
              z_list,
              objective_function,
              distribution_parameters)
g <- map2_dbl(x_list,
              z_list,
              constraint_function,
              distribution_parameters)

## data for plot with optimal points
plot_data <- bind_cols(x_grid, z_grid, tibble(f = f), tibble(g = g)) %>%
  mutate(z1 = factor(z1))
levels(plot_data$z1) <- c(`0` = TeX("$z_1 = 0$"), `1` = TeX("$z_1 = 1$"))

point_data <- tibble(z1 = factor(c(0,1)),
                     x1 = c(0.25, 0.75),
                     x2 = c(0.75, 0.25))
levels(point_data$z1) <- c(`0` = TeX("$z_1 = 0$"), `1` = TeX("$z_1 = 1$"))

## This is plot of data generating function
true_dgm_plot <- ggplot(plot_data, aes(x = x1, y = x2, z = f)) +
  facet_wrap(z1~., labeller = label_parsed) +
  geom_contour_filled(bins = 6) + 
  geomtextpath::geom_textcontour(aes(x = x1, y = x2, z = g),
               breaks = c(0.2, 0.5, 0.8, 1.2, 1.5),
               size = 2.5,
               color = "white",
               linetype = "dashed",
               #label = c("0.2", "0.5", "0.8", "1.2", "1.5"),
               inherit.aes = FALSE) +
  geom_point(data = point_data,
             aes(x = x1, y = x2), shape = 8, color = "white", inherit.aes = FALSE) +
  coord_cartesian(xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1)) +
  guides(fill = guide_legend(title = "", ncol = 1, byrow = FALSE)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) + #match scale of other plots
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  scale_fill_viridis_d(option = "D") +
  ylab(parse(text = TeX('$d_2$')))+
  xlab(parse(text = TeX('$d_1$')))+
  theme(panel.background = element_rect(fill='darkgrey'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        aspect.ratio = 1)

true_dgm_plot

#---------------------------------------------------------------------------#
set_null_device(cairo_ps) # this is required for the latex
final_plot <- plot_grid(true_dgm_plot,
                        plot_toxic_doses,
                        plot_exp_dose_units, 
                        plot_avg_rpsel,
                        nrow = 2,
                        byrow = TRUE,
                        labels = LETTERS[1:4],
                        scale = 1,
                        rel_heights = c(1,1.75))

ggsave2(filename = "bimodal_resp_het.eps",
        device = cairo_ps,
        plot = final_plot,
        width = 11,
        height = 8.5,
        units = "in",
        dpi = 800)









