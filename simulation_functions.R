## REQUIRES: Objective and toxicity functions to be defined w.r.t X's on UNIT CUBE! 
##  - the code does NOT rescale X's, so they are assumed to be on unit cube

pacman::p_load(hetGP, ggdist, dplyr, tidyr, purrr, foreach, doParallel, spacefillr, cli)

#------------- Below remains constant for each sim scenario -----------------------------#

# returns tibble of relevant performance metrics at dose combo x for strata z
quantify_performance <- function(z, #list
                                 best_x, #list
                                 mean_estimate_f,
                                 sd2_estimate_f,
                                 mean_estimate_g,
                                 sd2_estimate_g,
                                 true_f_at_best_x,
                                 true_g_at_best_x,
                                 distribution_parameters,
                                 x_grid_granularity,
                                 seed){
  
  # return missing values if no best x since implies constraints 
  # aren't met
  if(is.na(true_f_at_best_x) | is.na(true_g_at_best_x)){
    .res <- tibble(rpsel_true_f_at_optimal = NA_real_,
           rpsel_true_f_at_best_x = NA_real_,
           rpsel_true_g_at_optimal = NA_real_,
           rpsel_true_g_at_best_x = NA_real_,
           euclidean_distance_from_optimal = NA_real_,
           dose_units_from_optimal = NA_real_)
    return(.res)
  }

  # Find z from distribution_parameters which matches z of interest
  .index_z <- distribution_parameters$z %>%
    map(~ all(length(.) == length(flatten(z))) && all(. == flatten(z))) %>%
    list_c() %>%
    which()

  # Abort if z vector isn't part of data generating mechanism
  if(length(.index_z) == 0L){
    cli_abort(c(message = "Can't find z vector in {.var distribution_parameters$z}. Make sure it's the right length and has valid strata values.",
                "i" = "z vector: {z}",
                "i" = "x vector: {x}"))
  }
  
  .true_f_at_optimal <- distribution_parameters$true_f_at_optimal[[.index_z]]
  .true_g_at_optimal <- distribution_parameters$true_g_at_optimal[[.index_z]]
  .true_optimal_dose <- distribution_parameters$true_optimal_dose[[.index_z]]
  
  # posterior samples of f and g
  set.seed(seed)
  .posterior_samples_f <- rnorm(n = 10000,
                                mean = mean_estimate_f,
                                sd = sqrt(sd2_estimate_f))
  
  .posterior_samples_g <- rnorm(n = 10000,
                                mean = mean_estimate_g,
                                sd = sqrt(sd2_estimate_g))

  .rpsel_true_f_at_optimal <- sqrt(mean((.posterior_samples_f - .true_f_at_optimal)^2))
  .rpsel_true_f_at_best_x <- sqrt(mean((.posterior_samples_f - true_f_at_best_x)^2))
  
  .rpsel_true_g_at_optimal <- sqrt(mean((.posterior_samples_g - .true_g_at_optimal)^2))
  .rpsel_true_g_at_best_x <- sqrt(mean((.posterior_samples_g - true_g_at_best_x)^2))

  if(length(flatten(best_x)) != length(.true_optimal_dose)){
    cli_abort(c(message = "Dimension of {.var best_x} does not match dimension of {.var .true_optimal_dose}. ",
                "i" = "{.var best_x}: {best_x}",
                "i" = "{.var .true_optimal_dose}: {(.true_optimal_dose)}"))
  }
  
  .euclidean_distance_from_optimal <- map2(flatten(best_x), 
                                           .true_optimal_dose, 
                                           function(x,y) (x - y)^2) %>% 
                                      reduce(sum) %>% sqrt()

  .dose_units_from_optimal <- .euclidean_distance_from_optimal / x_grid_granularity

  tibble(rpsel_true_f_at_optimal = .rpsel_true_f_at_optimal,
         rpsel_true_f_at_best_x = .rpsel_true_f_at_best_x,
         rpsel_true_g_at_optimal = .rpsel_true_g_at_optimal,
         rpsel_true_g_at_best_x = .rpsel_true_g_at_best_x,
         euclidean_distance_from_optimal = .euclidean_distance_from_optimal,
         dose_units_from_optimal = .dose_units_from_optimal)
  
}

get_initial_data <- function(objective_function, 
                             constraint_function,
                             distribution_parameters,
                             true_noise_sd_f, 
                             true_noise_sd_g,
                             initial_n_doses,
                             x_grid_granularity,
                             seed,
                             reps_each_x,
                             x_generation_type = "dose_escalation",
                             z_enrollment_type = "balanced"){
  
  
  
  # There is a suggestion that for sequential tasks like BayesOpt, randomly assigning these initial
  # points can be better than space-filling
  # See: Distance-Distributed Design for Gaussian Process Surrogates
  # https://www.tandfonline.com/doi/full/10.1080/00401706.2019.1677269

  # Dimension must match mean vectors in data generating mechanism
  .dimension_x <- length(distribution_parameters$true_optimal_dose[[1]])
  .names_x <- paste0("x", 1:.dimension_x)

  if(x_generation_type == "dose_escalation"){
    .initial_x <- data.frame(matrix(rep(0, initial_n_doses * .dimension_x), ncol = .dimension_x, byrow = TRUE)) 
    names(.initial_x) <- .names_x
    
  } else if(x_generation_type == "random"){
    
    set.seed(seed)
    .lower <- c(0, 0)
    .upper <- c(1, 1)
    .initial_x <- data.frame(matrix(runif(initial_n_doses * .dimension_x, .lower, .upper),
                                    ncol = .dimension_x, byrow = TRUE)) 
    names(.initial_x) <- .names_x
    
  } else if (x_generation_type == "sobol") {
    
    .initial_x <- data.frame(generate_sobol_set(n = initial_n_doses, dim = .dimension_x, seed = seed))
    names(.initial_x) <- .names_x
  } 
  
  
  .initial_x <- .initial_x %>%
    map_dfr(plyr::round_any, accuracy = x_grid_granularity) %>%
    uncount(reps_each_x) %>%
    mutate(iteration = 0)
  
  .x_list <- .initial_x %>%
    select(-iteration) %>%
    transpose() %>%
    map(~ unlist(.), use.names = TRUE)
  
  if(z_enrollment_type == "balanced"){
  
    set.seed(seed)
    .z_row <- distribution_parameters$z %>%
      map(tibble::as_tibble_row) %>%
      list_rbind() %>%
      slice_sample(n = length(distribution_parameters$z), replace = FALSE) %>%
      uncount(reps_each_x)
    
    # If balance is selected, function forces at least one observation per stratum
    # which means initial_n_doses must be multiple of number of strata
    # e.g., 3 binary covaraites yields 8 strata, so 8, 16, 24, etc. initial doses
    # can be tried, which implies 1, 2, 3, etc. unique doses per strata
    
    .n_strata <- length(distribution_parameters$z)
    if(!(initial_n_doses %% .n_strata == 0)){
      cli_abort(c(message = "The number of initial doses {.var initial_n_doses} must be a multiple of the number of strata (i.e., length of {.var distribution_parameters$z}).",
                  "i" = "{.var initial_n_doses}: {initial_n_doses}",
                  "i" = "length of {.var distribution_parameters$z}: {length(distribution_parameters$z)}"))
    }
    
    .z_row <- bind_rows(replicate(initial_n_doses/.n_strata, 
                                  .z_row, simplify = FALSE))

    .z_list <- .z_row %>%
      transpose() %>%
      map(~ unlist(.), use.names = TRUE)

    .initial_x <- .initial_x %>%
      bind_cols(.z_row)


  } else if (z_enrollment_type == "random"){
    
    set.seed(seed)
    .z_row <- distribution_parameters$z %>%
      map(tibble::as_tibble_row) %>%
      list_rbind() %>%
      slice_sample(n = nrow(.initial_x), replace = TRUE)
    
    .z_list <- .z_row %>%
      transpose() %>%
      map(~ unlist(.), use.names = TRUE)
    
    .initial_x <- .initial_x %>%
      bind_cols(.z_row)
  }

  ## objective function
  .initial_f <-  map2_dbl(.x_list,
                      .z_list,
                      objective_function,
                      distribution_parameters)


  .initial_yf <- rnorm(n = length(.initial_f),
                      mean = .initial_f,
                      sd = true_noise_sd_f)
  
  ## constraint function
  .initial_g <-  map2_dbl(.x_list,
                          .z_list,
                          constraint_function,
                          distribution_parameters)
  
  
  .initial_yg <- rnorm(n = length(.initial_g),
                       mean = .initial_g,
                       sd = true_noise_sd_g)
  
  

  bind_cols(.initial_x,
            f = .initial_f,
            yf = .initial_yf,
            g = .initial_g,
            yg = .initial_yg)

}

get_data <- function(objective_function, 
                     constraint_function, 
                     distribution_parameters,
                     algorithm,
                     next_x,
                     evaluated_data,
                     iteration,
                     true_noise_sd_f, 
                     true_noise_sd_g, 
                     seed,
                     reps_each_x,
                     z_enrollment_type){
  
  #ensures different patterns of z over iterations within and across simulations
  set.seed(seed + iteration) 
  
  .strata <- distribution_parameters$z %>%
    map(tibble::as_tibble_row) %>%
    list_rbind()
  
  if(algorithm == "standard"){
    .x <- next_x %>%
      uncount(reps_each_x) 
    
    if (z_enrollment_type == "random"){
      #print(paste0("iteration: ", iteration))
      .x <- .x %>%
        bind_cols(slice_sample(.strata, n = nrow(.x), replace = TRUE))
      
    } else if (z_enrollment_type == "balanced"){
      .order_z <- evaluated_data$z %>%
        map(tibble::as_tibble_row) %>%
        list_rbind() %>%
        distinct() %>% #note that this preserves the order of z in evaluated_data
        mutate(order = row_number())

      .last_index <- nrow(evaluated_data)
      .last_val_z <- evaluated_data$z[.last_index] %>% map(tibble::as_tibble_row) %>% list_rbind()
      .last_val_z_order <- .order_z %>%
        right_join(.last_val_z, by = names(.last_val_z))
      
      if(nrow(.x) == 1){
        .new_row_order <- if_else(.last_val_z_order$order == max(.order_z$order),
                                  min(.order_z$order), 
                                  .last_val_z_order$order + 1L) %>% tibble::as_tibble_col(column_name = "order")
      } else {
        
        .new_row_order <- c(.last_val_z_order$order, rep(NA_integer_, nrow(.x)))
        for(i in 1:(nrow(.x))){
          .new_row_order[i+1] <-  if_else(.new_row_order[i] == max(.order_z$order),
                                          min(.order_z$order), 
                                          .new_row_order[i] + 1L)
        }
        .new_row_order <- .new_row_order[-1] %>% tibble::as_tibble_col(column_name = "order")
      }

      
      .new_z <- .new_row_order %>% left_join(.order_z, by = "order") %>% select(-order)

      .x <- .x %>%
        bind_cols(.new_z)

    }
    
  } else if (algorithm == "personalized"){
    
    if (z_enrollment_type == "random"){
      # randomly choose the next z status, single dose only
      .x <- slice_sample(.strata, n = 1, replace = TRUE) %>%
        left_join(next_x, by = names(.strata)) %>% # these won't match necessarily
        uncount(reps_each_x)
      
    } else if (z_enrollment_type == "balanced"){
      .order_z <- evaluated_data$z %>%
        map(tibble::as_tibble_row) %>%
        list_rbind() %>%
        distinct()
      .x <- .order_z %>%
        left_join(next_x, by = names(.order_z)) %>%
        uncount(reps_each_x)
    }
    
  }

  .x_active_strata <- .x %>% filter(af_early_stop == 0, constraint_early_stop == 0) 
  if(nrow(.x_active_strata) == 0){
    .x <- .x %>%
      mutate(z = .x %>%
               select(starts_with("z")) %>%
               transpose() %>%
               map(~ unlist(.), use.names = TRUE),
             f = NA_real_,
             yf = NA_real_,
             g = NA_real_,
             yg = NA_real_)
    return(.x)
  }

  foreach(i = 1:nrow(.x)) %do% {
    .x_list <- .x %>%
      select(starts_with("x")) %>%
      transpose() %>%
      map(~ unlist(.), use.names = TRUE)
    .z_list <- .x %>%
      select(starts_with("z")) %>%
      transpose() %>%
      map(~ unlist(.), use.names = TRUE)
  }
  
  # objective function
  .f <-  map2_dbl(.x_list,
                  .z_list,
                  objective_function,
                  distribution_parameters)
  
  # generate NAs if .f is NA, as intended
  .yf <- suppressWarnings(rnorm(n = length(.f),
              mean = .f,
              sd = true_noise_sd_f))
  
  # constraint function
  .g <-  map2_dbl(.x_list,
                  .z_list,
                  constraint_function,
                  distribution_parameters)
  
  # generate NAs if .g is NA, as intended
  .yg <- suppressWarnings(rnorm(n = length(.g),
              mean = .g,
              sd = true_noise_sd_g))
  
  bind_cols(.x,
            f = .f,
            yf = .yf,
            g = .g,
            yg = .yg,
            tibble(z = .z_list)) %>%
    mutate(f = if_else(af_early_stop == 1 | constraint_early_stop == 1, NA_real_, f),
           yf = if_else(af_early_stop == 1 | constraint_early_stop == 1, NA_real_, yf),
           g = if_else(af_early_stop == 1 | constraint_early_stop == 1, NA_real_, g),
           yg = if_else(af_early_stop == 1 | constraint_early_stop == 1, NA_real_, yg))
  
}

#----------------- Acquisition Functions -------------------------#

# Must return VECTOR of values of acquisition function corresponding to each row of 
# optimization grid .new_x_grid

# Constrained Expected Improvement
  # crit_EI does NOT use the noise variance/nugget when performing EI calculation
  # since original EI was proposed in deterministic setting
  # It also doesn't rescale to original output scale so make sure
  # cst uses something on scale which rescale_data() produces

CEI <- function(.predictions_at_x, 
                .observed_mean_yg,
                .observed_sd_yg,
                .best_observed_f, #on STANDARDIZED SCALE
                .new_x_grid, 
                .gp_mod_f){
  
  # This does not rescale f, so best_observed_f should be on STANDARDIZED scale

  .ei <- crit_EI(x = as.matrix(.new_x_grid),
          model = .gp_mod_f,
          cst = .best_observed_f)
  .constraint_threshold_g <- unique(.predictions_at_x$constraint_threshold_g)
  # the constraint threshold needs to be standardized using the values from
  # rescale_data() for g 
  .rescaled_constraint <- (.constraint_threshold_g - .observed_mean_yg) / .observed_sd_yg
  
  .ei * pnorm((.rescaled_constraint - .predictions_at_x$mean_g)/sqrt(.predictions_at_x$sd2_g))

}


  #---------- Stratum-specific early stopping ----------------#

determine_stratum_stop <- function(next_x, 
                                   evaluated_data, 
                                   af_stop_threshold, 
                                   iteration){
  
  .lag_1_index <- max(unique(evaluated_data$iteration))
  .lag_2_index <- .lag_1_index - 1
  .lag_3_index <- .lag_2_index - 1

  # AF check goes back only to two places since check in advance of third
    if(next_x$max_af_value < af_stop_threshold){
      .below_af_stop_threshold <- 1
    } else {
      .below_af_stop_threshold <- 0
    }
    
    .below_af_stop_threshold_lag_1 <- evaluated_data %>% 
      filter(iteration == .lag_1_index) %>%
      pull(below_af_stop_threshold) %>% unique()
    .below_af_stop_threshold_lag_2 <- evaluated_data %>% 
      filter(iteration == .lag_2_index) %>%
      pull(below_af_stop_threshold) %>% unique()
    
    .sum_below_af_stop_threshold <- sum(.below_af_stop_threshold,
                                        .below_af_stop_threshold_lag_1,
                                        .below_af_stop_threshold_lag_2)
    
    # feasibility goes back three places since relies on feasibility at particular iteration
    .no_feasible_x_in_stratum_lag_1 <- evaluated_data %>% 
      filter(iteration == .lag_1_index) %>%
      pull(no_feasible_x_in_stratum) %>% unique()
    .no_feasible_x_in_stratum_lag_2 <- evaluated_data %>% 
      filter(iteration == .lag_2_index) %>%
      pull(no_feasible_x_in_stratum) %>% unique()
    .no_feasible_x_in_stratum_lag_3 <- evaluated_data %>% 
      filter(iteration == .lag_3_index) %>%
      pull(no_feasible_x_in_stratum) %>% unique()
    
    .sum_no_feasible_x_in_stratum <- sum(.no_feasible_x_in_stratum_lag_1,
                                         .no_feasible_x_in_stratum_lag_2,
                                         .no_feasible_x_in_stratum_lag_3)
    
    .previous_af_early_stop <- evaluated_data %>% 
      filter(iteration == .lag_1_index) %>%
      pull(af_early_stop) %>%
      unique()
    
    .previous_constraint_early_stop <- evaluated_data %>% 
      filter(iteration == .lag_1_index) %>%
      pull(constraint_early_stop) %>%
      unique()
    
    bind_cols(next_x,
              tibble(iteration = iteration,
                     below_af_stop_threshold = .below_af_stop_threshold,
                     af_early_stop = if_else(.previous_af_early_stop == 1 |
                                               .sum_below_af_stop_threshold == 3, 1, 0),
                     constraint_early_stop = if_else(.previous_constraint_early_stop == 1 |
                                                     .sum_no_feasible_x_in_stratum == 3, 1, 0)))
    
}
  
## Standard dose-finding simulation function 
dose_finding <- function(simulation,
                        algorithm,
                        objective_function, 
                        constraint_function,
                        constraint_threshold_g, #should take form of tibble, no z if standard, z if personalized
                        distribution_parameters,
                        true_noise_sd_f, 
                        true_noise_sd_g,
                        initial_n_doses, 
                        x_grid_granularity, 
                        seed, 
                        reps_each_x,
                        max_sample_size,
                        acquisition_function,
                        x_generation_type,
                        z_enrollment_type,
                        af_stop_threshold,
                        constraint_probability_threshold,
                        rho,
                        debug = FALSE) {
  
  # Note that conditions are introduced throughout to separate the different
  # functional behavior between the standard and personalized algorithms while 
  # preserving their shared functional behavior
  
  #----------------------------------------------------------------------------#
  #------------------ Initial Step --------------------------------------------#
  #----------------------------------------------------------------------------#
  
  
  #---------------- Get initial data ------------------------#
  .initial_data <- get_initial_data(objective_function, 
                                    constraint_function,
                                    distribution_parameters,
                                    true_noise_sd_f, 
                                    true_noise_sd_g,
                                    initial_n_doses, 
                                    x_grid_granularity, 
                                    seed, 
                                    reps_each_x,
                                    x_generation_type,
                                    z_enrollment_type)

  .initial_data_names <- names(.initial_data)
  .x_columns_names <- .initial_data_names[stringr::str_detect(.initial_data_names, pattern = "x\\d+")]
  
  if(algorithm == "standard"){

    .z_columns_names <- NULL
    .initial_matrix_x <- as.matrix(.initial_data[, .x_columns_names])
  
  } else if (algorithm == "personalized"){

    .z_columns_names <- .initial_data_names[stringr::str_detect(.initial_data_names, pattern = "z\\d+")] 
    .initial_matrix_x <- as.matrix(.initial_data[, c(.x_columns_names, .z_columns_names)])
 
  }

  # objective function
  .initial_vector_f <- pull(.initial_data[,c("f")]) 
  .initial_vector_yf <- pull(.initial_data[,c("yf")])
  
  # constraint function
  .initial_vector_g <- pull(.initial_data[,c("f")]) 
  .initial_vector_yg <- pull(.initial_data[,c("yg")])
  
  #--------------- Get initial fit of GPs ---------------------#
  # Standardize response, inputs are already on unit cube (THIS IS ASSUMED) so 
  # do not need to be rescaled
  
  # objective function
  .rescaled_initial_data_f <- find_reps(X = .initial_matrix_x,
                                     Z = .initial_vector_yf,
                                     rescale = FALSE,
                                     normalize = TRUE)

  .observed_mean_yf <- .rescaled_initial_data_f$outputStats[1] 
  .observed_sd_yf <- sqrt(.rescaled_initial_data_f$outputStats[2])
  

  # constraint function
  .rescaled_initial_data_g <- find_reps(X = .initial_matrix_x,
                                        Z = .initial_vector_yg,
                                        rescale = FALSE,
                                        normalize = TRUE)
  
  .observed_mean_yg <- .rescaled_initial_data_g$outputStats[1] 
  .observed_sd_yg <- sqrt(.rescaled_initial_data_g$outputStats[2])
  
  #objective function
  # Note for mleHomGP that Z is the standardized response, NOT additional covariates 
  
  # Only single point under dose_escalation scheme, so set ML estimates of hyperparameters 
  # to fixed value until enough data to converge 
  if(x_generation_type == "dose_escalation"){
    
    .gp_mod_f <- suppressWarnings(mleHomGP(.rescaled_initial_data_f,
                           .rescaled_initial_data_f$Z,
                           lower = rep(sqrt(.Machine$double.eps), ncol(.initial_matrix_x)),
                           upper = rep(sqrt(ncol(.initial_matrix_x)), ncol(.initial_matrix_x)),
                           known = list(g = matrix(.observed_sd_yf^2), 
                                        theta = rep(sqrt(ncol(.initial_matrix_x))/2, ncol(.initial_matrix_x)))
                           ))
    
  } else {
    .gp_mod_f <-  mleHomGP(.rescaled_initial_data_f,
                           .rescaled_initial_data_f$Z,
                           settings = list(pgtol = 1e-05))
  }
  

              # Note the following: 
              #
              # 1) Restriction of length-scales theta to (0, sqrt(x dimension)) comes from
              # Distance-Distributed Design for Gaussian Process Surrogates
              # https://www.tandfonline.com/doi/full/10.1080/00401706.2019.1677269
              # Related - data on unit cube won't be very informative outside this range
              #
              # 2) pgtol = 1e-05, for convergence of optimizer as specified in L-BFGS-B documentation 
              #     https://github.com/jonathanschilling/L-BFGS-B/blob/gh-pages/L-BFGS-B.pdf

  
    # constraint function
  if(x_generation_type == "dose_escalation"){
    .gp_mod_g <-  suppressWarnings(mleHomGP(.rescaled_initial_data_g,
                                            .rescaled_initial_data_g$Z,
                                            lower = rep(sqrt(.Machine$double.eps), ncol(.initial_matrix_x)),
                                            upper = rep(sqrt(ncol(.initial_matrix_x)), ncol(.initial_matrix_x)),
                                            known = list(g = matrix(.observed_sd_yg^2), 
                                                         theta = rep(sqrt(ncol(.initial_matrix_x))/2, ncol(.initial_matrix_x)))))
  } else {
  .gp_mod_g <-  mleHomGP(.rescaled_initial_data_g,
                         .rescaled_initial_data_g$Z,
                         settings = list(pgtol = 1e-05))
  }
  
  # If initial model doesn't converge for random/sobol x selection, abort
  if(stringr::str_detect(.gp_mod_f$msg, "ERROR")){
    cli_abort(c(message = "Initial GP model for f didn't converge",
                "i" = "seed: {seed}",
                "i" = "true_noise_sd_f: {true_noise_sd_f}",
                "i" = "initial_n_doses: {initial_n_doses}",
                "i" = "reps_each_x: {reps_each_x}",
                "i" = "iteration: {unique(.initial_data$iteration)}",
                "i" = "message: {(.gp_mod_f$msg)}"))
  }
  
  if(stringr::str_detect(.gp_mod_g$msg, "ERROR")){
    cli_abort(c(message = "Initial GP model for g didn't converge",
                "i" = "seed: {seed}",
                "i" = "true_noise_sd_g: {true_noise_sd_g}",
                "i" = "initial_n_doses: {initial_n_doses}",
                "i" = "reps_each_x: {reps_each_x}",
                "i" = "iteration: {unique(.initial_data$iteration)}",
                "i" = "message: {(.gp_mod_g$msg)}"))
  }
  
  
  
if(debug == TRUE){
  print("made it upt to grid_x")
}

  
  
  #------------- Grid for optimization -------------------#
  # Full grid for random/sobol x generation
  .grid_x <- .x_columns_names %>%
    map(~tibble(!!.x := seq(0,1, x_grid_granularity))) %>%
    list_cbind() %>%
    reduce(expand_grid) %>%
    set_names(.x_columns_names)
  
  
  
if(debug == TRUE){
print("made it before joining constraint threshold")  
}
  
  
  
  if(algorithm == "standard"){
    .new_x_grid <- .grid_x %>% bind_cols(constraint_threshold_g)
  } else if (algorithm == "personalized"){
    .grid_z <- distribution_parameters$z %>%
      map(tibble::as_tibble_row) %>%
      list_rbind()
    .new_x_grid <- expand_grid(.grid_x, .grid_z) %>%
      left_join(constraint_threshold_g, by = names(.grid_z))
  }
  
if(debug == TRUE){
  print("made it before dose escalation") 
}
  
  ### Under dose escalation, can only consider points contained in expanding 
  ### dose combination region
  ### - can generalize this to include different rho for strata, future work
  if(x_generation_type == "dose_escalation"){
    # take intersection of unit hypercube and
    # lower half-space for iteration 0
    .iteration <- 0
    .new_x_grid <- .new_x_grid %>%
      rowwise() %>%
      mutate(permitted_region = if_else(sum(c_across(starts_with("x"))) <= rho * .iteration, 1, 0)) %>%
      ungroup() %>%
      filter(permitted_region == 1) %>%
      select(-permitted_region)
    
  }
 
  
  
  #------------- Obtain best x and f ---------------------------------------------#
  #------------- subject to satisfying constraint threshold on g -----------------#
  .predictions_f_at_x <- predict(.gp_mod_f,
                            x = as.matrix(.new_x_grid %>% select(-constraint_threshold_g))) 
  .predictions_g_at_x <- predict(.gp_mod_g,
                                 x = as.matrix(.new_x_grid %>% select(-constraint_threshold_g))) 
  
  .z_columns_names_syms <-  if(is.null(.z_columns_names)) 'DROP_COLUMN' else syms(.z_columns_names)

  
  
if(debug == TRUE){
  print("made it up to minimizer")
}

  
  
  .minimizer <- .new_x_grid %>%
    # convert back to original scale
    mutate(posterior_mean_f_at_best_x = .predictions_f_at_x$mean*.observed_sd_yf + .observed_mean_yf,
           posterior_sd2_f_at_best_x = .predictions_f_at_x$sd2*(.observed_sd_yf)^2, 
           posterior_mean_g_at_best_x = .predictions_g_at_x$mean*.observed_sd_yg + .observed_mean_yg,
           posterior_sd2_g_at_best_x = .predictions_g_at_x$sd2*(.observed_sd_yg)^2) %>%
    mutate(probability_satisfying_constraint_g = pmap_dbl(list(.$posterior_mean_g_at_best_x,
                                                               sqrt(.$posterior_sd2_g_at_best_x),
                                                               .$constraint_threshold_g),
                                                      function(x,y,z) pnorm(q = z,
                                                                          mean = x,
                                                                          sd = y)),
           feasible_best_x = if_else(probability_satisfying_constraint_g > 
                                  constraint_probability_threshold, 1, 0)) %>%
    group_by(!!!.z_columns_names_syms)
  
  .feasible_x_check <- .minimizer %>%
      summarize(no_feasible_x_in_stratum = if_else(all(feasible_best_x == 0), 1L, 0L)) 

  ## Three situations 
  ## 1) all infeasible - then select point where prob_satisfying_constraint is highest (Gelbart et al 2014)
  ##     need to ensure enough data devoted to this initial step otherwise run the chance of 
  ##     terminating too early
  ## 2) all feasible - filter to feasible points and select minimizer of posterior mean of objective
  ## 3) mixed - do both above for strata with infeasible/feasible points respectively
  
      .minimizer_temp <- .minimizer %>%
        left_join(.feasible_x_check, by = .feasible_x_check %>% select(-no_feasible_x_in_stratum) %>% names())
      
      if(all(.feasible_x_check$no_feasible_x_in_stratum == 1L)){
        .minimizer <- .minimizer_temp %>%
          filter(probability_satisfying_constraint_g == max(probability_satisfying_constraint_g)) %>%
          slice_sample(n = 1) 
        
      } else if(all(.feasible_x_check$no_feasible_x_in_stratum == 0L)){
        .minimizer <- .minimizer_temp %>%
          filter(feasible_best_x == 1) %>%
          filter(posterior_mean_f_at_best_x == min(posterior_mean_f_at_best_x)) %>%
          slice_sample(n = 1)
        
      } else {
        .minimizer_feasible <-  .minimizer_temp %>%
          filter(no_feasible_x_in_stratum == 0) %>%
          filter(feasible_best_x == 1) %>%
          filter(posterior_mean_f_at_best_x == min(posterior_mean_f_at_best_x)) %>%
          slice_sample(n = 1)
        
        .minimizer_no_feasible <-  .minimizer_temp %>%
          filter(no_feasible_x_in_stratum == 1) %>%
          filter(probability_satisfying_constraint_g == max(probability_satisfying_constraint_g)) %>%
          slice_sample(n = 1)
        
        .minimizer <- .minimizer_feasible %>%
          bind_rows(.minimizer_no_feasible) #%>%
        
      }
    
      
      
if(debug == TRUE){
  print("made it past three part section") 
}  
      
      
     
      .minimizer <- .minimizer %>%
        ungroup() %>%
        select(-any_of("\"DROP_COLUMN\"")) %>%
        rename_with(~ paste0("best_", .x_columns_names), starts_with("x"))
      .x_columns <- .minimizer %>% select(starts_with("best_x"))
      .minimizer <- .minimizer %>%
        mutate(best_x = .x_columns %>% transpose() %>% map(~ unlist(.), use.names = TRUE))
          
      if(algorithm == "standard"){
        # get z columns from initial data since not included under standard algorithm
        .z_columns <- .initial_data %>%
          select(!!!syms(.initial_data_names[stringr::str_detect(.initial_data_names, pattern = "z\\d+")])) %>%
          distinct()
        .minimizer <- .minimizer %>%
          bind_cols(.z_columns)
      }
    
      .z_columns <- .minimizer %>% select(starts_with("z")) 
      .minimizer <- .minimizer %>%
        mutate(z = .z_columns %>% transpose() %>% map(~ unlist(.), use.names = TRUE),
               max_af_value = NA_real_,
               gp_message_f = .gp_mod_f$msg,
               gp_message_g = .gp_mod_g$msg) %>%
        mutate(true_f_at_best_x = map2_dbl(.$best_x, .$z, objective_function, distribution_parameters),
               true_g_at_best_x = map2_dbl(.$best_x, .$z, constraint_function, distribution_parameters))
      .evaluated_data <- .initial_data %>%
        left_join(.minimizer, by = names(.z_columns))
      
      #----------- Quantify performance metrics ---------------#
      # don't save posterior samples as memory requirements are 
      # too large given number of iterations/simulations
      .performance <- foreach(i = 1:nrow(.evaluated_data), .combine = 'rbind') %do% {
        quantify_performance(z = .evaluated_data$z[i],
                            best_x = .evaluated_data$best_x[i],
                            mean_estimate_f = .evaluated_data$posterior_mean_f_at_best_x[i],
                            sd2_estimate_f = .evaluated_data$posterior_sd2_f_at_best_x[i],
                            mean_estimate_g = .evaluated_data$posterior_mean_g_at_best_x[i],
                            sd2_estimate_g = .evaluated_data$posterior_sd2_g_at_best_x[i],
                            true_f_at_best_x = .evaluated_data$true_f_at_best_x[i],
                            true_g_at_best_x = .evaluated_data$true_g_at_best_x[i],
                            distribution_parameters = distribution_parameters,
                            x_grid_granularity = x_grid_granularity,
                            seed = seed)
      }

      .evaluated_data <- bind_cols(.evaluated_data,
                                   .performance,
                                   af_early_stop = 0, #at specific iteration
                                   below_af_stop_threshold = 0,
                                   constraint_early_stop = 0) # at specific iteration
      
      
if(debug == TRUE){
  print("made it past initial evaluated data") 
}

      

  #----------------------------------------------------------------------#
  #------------------ Dose Finding --------------------------------------#
  #----------------------------------------------------------------------#
  
  
  
  #----------- Sanity check ------------------------#
  # Determine max number of algorithm iterations and ensure it 
  # is an integer and non-negative
  if(algorithm == "personalized" & z_enrollment_type == "balanced"){
    .max_n_iterations <- (max_sample_size - initial_n_doses*reps_each_x)/
      (reps_each_x*length(distribution_parameters$z))
  } else {
    .max_n_iterations <- (max_sample_size - initial_n_doses*reps_each_x)/reps_each_x
  }
  
  .integer_check <- floor(.max_n_iterations) - .max_n_iterations == 0
  .non_negativity_check <- .max_n_iterations >= 0

      if(.integer_check == FALSE | .non_negativity_check == FALSE){
        cli_abort(c(message = "For dose finding, max_n_iterations is not an integer or is negative",
                    "i" = "For personalized algortihm with balanced z_enrollment_type, max_n_iterations = (max_sample_size - initial_n_doses*reps_each_x)/(reps_each_x*length(distribution_parameters$z))",
                    "i" = "For everything else, max_n_iterations = (max_sample_size - initial_n_doses*reps_each_x)/reps_each_x",
                    "i" = "max_n_iterations = {(.max_n_iterations)}"))
      }
  
  
  
  #----------- Initialize early stopping conditions -----------#
  
  .n_enrolled <- nrow(.evaluated_data %>% filter(!is.na(f))) 
  .all_strata_stopped <- .evaluated_data %>% 
    group_by(!!!.z_columns_names_syms) %>%
    slice(n()) %>% ungroup() %>% select(af_early_stop, constraint_early_stop) %>% mutate(any_stop = if_else(af_early_stop == 1 | constraint_early_stop == 1, 1, 0)) #no_feasible_x_in_stratum
  .all_strata_stopped_logical <- all(.all_strata_stopped$any_stop == 1) 
  .iteration <- 1

  #----------- Begin search algorithm -------------------------#
  while(.n_enrolled < max_sample_size & .all_strata_stopped_logical == FALSE){

    if(.max_n_iterations == 0) {
      break
    }

    
    
if(debug == TRUE){
  print(paste0("begin iteration: ", .iteration))  
}    
 
 
    #---- Determine next dose to evaluate --------------------------#
    # For personalized setting, calculated for each stratum and then later 
    # filtered if needed for z_generation_type = "random"
    # recall .z_columns_names is NULL if standard algorithm is used
    # posterior_mean_f must be RE-STANDARDIZED to be used in crit_EI function
   
     .best_observed_f <- .minimizer %>%
      select(any_of(c("posterior_mean_f_at_best_x", .z_columns_names))) %>%
      distinct() %>%
      mutate(posterior_mean_f_at_best_x = (posterior_mean_f_at_best_x - .observed_mean_yf) / .observed_sd_yf)
     
     ### Update .new_x_grid and obtain predictions if 
     ### x_generation_type is dose_escalation
     if(x_generation_type == "dose_escalation"){
       
       if(algorithm == "standard"){
         .new_x_grid <- .grid_x %>% bind_cols(constraint_threshold_g)
       } else if (algorithm == "personalized"){
         .grid_z <- distribution_parameters$z %>%
           map(tibble::as_tibble_row) %>%
           list_rbind()
         .new_x_grid <- expand_grid(.grid_x, .grid_z) %>%
           left_join(constraint_threshold_g, by = names(.grid_z))
       }
       
     # take intersection of unit hypercube and
     # lower half-space for iteration
         .new_x_grid <- .new_x_grid %>%
           rowwise() %>%
           mutate(permitted_region = if_else(sum(c_across(starts_with("x"))) <= rho * .iteration, 1, 0)) %>%
           ungroup() %>%
           filter(permitted_region == 1) %>%
           select(-permitted_region)
         
      ## NOW REMOVE PREVIOUSLY EXPLORED POINTS FOR FIRST ITERATIONS
      ## Until D is [0,1]^J
         if(rho * .iteration <= length(.x_columns_names)){
           #previously tried points 
           
           if(algorithm == "personalized"){
             .already_evaluated_x <- .evaluated_data %>%
               select(any_of(c(num_range("x", 1:20), num_range("z", 1:20))), constraint_threshold_g) %>%
               distinct()
           } else {
             .already_evaluated_x <- .evaluated_data %>%
               select(num_range("x", 1:20), constraint_threshold_g) %>%
               distinct()
           }
           
      
           .new_x_grid <- setdiff(.new_x_grid,  .already_evaluated_x)
             
         }
         
         .predictions_f_at_x <- predict(.gp_mod_f,
                                        x = as.matrix(.new_x_grid %>% select(-constraint_threshold_g))) 
         .predictions_g_at_x <- predict(.gp_mod_g,
                                        x = as.matrix(.new_x_grid %>% select(-constraint_threshold_g))) 
     }

     
     
if(debug == TRUE){
  print("made it up to af values")   
}
  
     
     
    # return tibble with list columns for consistency across algorithms and AFs    
    .af_values <- if(acquisition_function == "CEI") { # Constrained Expected Improvement

      # can join the stratum-specific g from distribution_parameters object here
        .predictions_af <- bind_cols(.new_x_grid,
                                  mean_f = .predictions_f_at_x$mean,
                                  sd2_f = .predictions_f_at_x$sd2,
                                  nugs_f = .predictions_f_at_x$nugs,
                                  mean_g = .predictions_g_at_x$mean,
                                  sd2_g = .predictions_g_at_x$sd2,
                                  nugs_g = .predictions_g_at_x$nugs) 
        
          if(algorithm == "standard"){
            
            
            list(CEI(.predictions_af,
                     .observed_mean_yg,
                     .observed_sd_yg,
                     .best_observed_f$posterior_mean_f_at_best_x, # this is standardized above
                     .new_x_grid %>% select(-constraint_threshold_g), 
                     .gp_mod_f))
            
          } else if (algorithm == "personalized"){
            
            
            
            foreach(i = 1:nrow(.best_observed_f)) %do%{
              
              .best_observed_f_stratum <- .best_observed_f %>%
                filter(row_number() == i)

              CEI(.predictions_af %>%
                    right_join(.best_observed_f_stratum %>%
                                 select(-posterior_mean_f_at_best_x), by = .z_columns_names),
                  .observed_mean_yg,
                  .observed_sd_yg,
                  .best_observed_f_stratum %>%
                    pull(posterior_mean_f_at_best_x),
                  .new_x_grid %>%
                    select(-constraint_threshold_g) %>%
                    right_join(.best_observed_f_stratum %>%
                                 select(-posterior_mean_f_at_best_x), by = .z_columns_names), 
                  .gp_mod_f)
              
            }
          }
    }
    
    # returns NA if best_x is missing
    .af_values <- bind_cols(.best_observed_f %>% select(-posterior_mean_f_at_best_x),
              tibble(af_values = .af_values))
    
    
    
    #------------ Find maximizer of AF and determine next x ---------------------#

    .index_max_af <- map_dbl(.af_values %>% select(af_values) %>% flatten(), 
                             function(x) ifelse(all(is.na(x)), NA_real_, which.max(x)))
    # don't use dplyr::if_else() here since it will try to evaluate false condition to ensure
    # type consistency so in event of NA vector which.max will return integer(0) and throw error
    
    .max_af_value <- map2_dbl(.af_values %>% select(af_values) %>% flatten(),
                              .index_max_af,
                              function(x,y) if(is.na(y)) NA_real_ else x[y])
    
    
if(debug == TRUE){
  print("made it past max af value")   
}
    
    if(algorithm == "standard"){
      
      .next_x <- .new_x_grid[.index_max_af, ]
      
    } else if (algorithm == "personalized"){
      
      .next_x <- foreach(i = 1:nrow(.af_values), .combine = 'rbind') %do% {
        .new_x_grid_subset <- .new_x_grid %>%
          right_join(.af_values %>%
                       filter(row_number() == i) %>%
                       select(-af_values), by = .z_columns_names)
        
        if(is.na(.index_max_af[i])){ #make x NA if no feasible points in stratum
          .new_x_grid_subset[1,] %>%
            mutate(across(starts_with("x"), ~ if_else(is.double(.), NA_real_, NA_real_)))
          
        } else {
          .new_x_grid_subset[.index_max_af[i],]
        }
        
      }
      
    }
    
    .next_x <- bind_cols(.next_x,
                         tibble(max_af_value = .max_af_value))
 
    
    
if(debug == TRUE){
  print("made it past next_x") 
}   
   
    
    #------------ Determine whether or not to stop early by stratum -------------#
    # Stop if AF value falls below threshold (J+1) times 
    # Or stop if no feasible doses (J+1) times
    # J is dimension of dosing agents x
    # e.g., J = 2 dose agents x1 and x2
    # keep this fixed at J = 2 but generalize in future if needed
    # Note: once stratum is stopped, enrollment CANNOT start again
    # this behavior can be changed by modifying determine_stratum_stop()

    .next_x <- if(algorithm == "standard"){
      determine_stratum_stop(.next_x, 
                             .evaluated_data, 
                             af_stop_threshold, 
                             .iteration)
    
    } else if (algorithm == "personalized"){

          foreach(i = 1:nrow(.next_x), .combine = 'rbind') %do% {
            determine_stratum_stop(.next_x[i, ],
                                   .evaluated_data %>%
                                     right_join(.next_x[i, ] %>%
                                                  select(!!!.z_columns_names_syms),
                                                by = .z_columns_names),
                                   af_stop_threshold,
                                   .iteration)
        }
  
    }
    
    
    
if(debug == TRUE){
  print("made it past stopping rule determination")   
}
 
    
    #----------- Evaluate new doses ------------------------------#
    .new_data <- get_data(objective_function = objective_function, 
                          constraint_function = constraint_function,
                          distribution_parameters = distribution_parameters,
                          algorithm = algorithm,
                          next_x = .next_x,
                          evaluated_data = .evaluated_data,
                          iteration = .iteration,
                          true_noise_sd_f = true_noise_sd_f, 
                          true_noise_sd_g = true_noise_sd_g, 
                          seed = seed,
                          reps_each_x = reps_each_x,
                          z_enrollment_type = z_enrollment_type)

    
if(debug == TRUE){
  print("made it past new_data generation")  
}
  
    
    
    #----------- Fit new GP model ---------------------------------#
    # keep only rows where dose-finding is active - no af_early_stop and
    # no constraint_early_stop
    .model_data <- bind_rows(.evaluated_data,
                             .new_data) %>%
      filter(af_early_stop == 0, constraint_early_stop == 0) 
    
    
    if(algorithm == "standard"){

      .matrix_x <- as.matrix(.model_data[, .x_columns_names])

    } else if (algorithm == "personalized"){

      .matrix_x <- as.matrix(.model_data[, c(.x_columns_names, .z_columns_names)])

    }
    
    
    # objective
    .vector_f <- pull(.model_data[,c("f")])
    .vector_yf <- pull(.model_data[,c("yf")])
    
    # constraint
    .vector_g <- pull(.model_data[,c("g")])
    .vector_yg <- pull(.model_data[,c("yg")])
    
    
    .rescaled_model_data_f <- find_reps(X = .matrix_x,
                                        Z = .vector_yf,
                                        rescale = FALSE,
                                        normalize = TRUE)
    .observed_mean_yf <- .rescaled_model_data_f$outputStats[1] 
    .observed_sd_yf <- sqrt(.rescaled_model_data_f$outputStats[2])
    
    .rescaled_model_data_g <- find_reps(X = .matrix_x,
                                        Z = .vector_yg,
                                        rescale = FALSE,
                                        normalize = TRUE)
    .observed_mean_yg <- .rescaled_model_data_g$outputStats[1] 
    .observed_sd_yg <- sqrt(.rescaled_model_data_g$outputStats[2])

    
   
    # Update each model only if new one converges
    # Otherwise keep old
    # See section 4.2 of "Noisy kriging-based optimization methods: A unified implementation within the DiceOptim package"
    # https://www.sciencedirect.com/science/article/abs/pii/S0167947313001205
    
    # objective
    .gp_mod_temp_f <-  try(mleHomGP(.rescaled_model_data_f,
                              .rescaled_model_data_f$Z,
                              # lower = rep(sqrt(.Machine$double.eps), 2),
                              # upper = rep(sqrt(2), 2),
                              settings = list(pgtol = 1e-05)), silent = TRUE)


    # if error in fitting model, use old one
    if(class(.gp_mod_temp_f) != "homGP"){

      .gp_mod_f <- suppressWarnings(mleHomGP(.rescaled_model_data_f,
                                  .rescaled_model_data_f$Z,
                                  lower = rep(sqrt(.Machine$double.eps), ncol(.matrix_x)),
                                  upper = rep(sqrt(ncol(.matrix_x)), ncol(.matrix_x)),
                                  known = list(g = matrix(.observed_sd_yf^2),
                                               theta = .gp_mod_f$theta)))
      
    } else if (stringr::str_detect(.gp_mod_temp_f$msg, "ERROR")){
      .gp_mod_f <- suppressWarnings(mleHomGP(.rescaled_model_data_f,
                          .rescaled_model_data_f$Z,
                          lower = rep(sqrt(.Machine$double.eps), ncol(.matrix_x)),
                          upper = rep(sqrt(ncol(.matrix_x)), ncol(.matrix_x)),
                          known = list(g = .gp_mod_f$g,
                                       theta = .gp_mod_f$theta)))
    } else {
      .gp_mod_f <- .gp_mod_temp_f
    }
    

    # constraint surface
    .gp_mod_temp_g <-  try(mleHomGP(.rescaled_model_data_g,
                                .rescaled_model_data_g$Z,
                                settings = list(pgtol = 1e-05)), silent = TRUE)
    
    if(class(.gp_mod_temp_g) != "homGP"){
      .gp_mod_g <- suppressWarnings(mleHomGP(.rescaled_model_data_g,
                                .rescaled_model_data_g$Z,
                                lower = rep(sqrt(.Machine$double.eps), ncol(.matrix_x)),
                                upper = rep(sqrt(ncol(.matrix_x)), ncol(.matrix_x)),
                                known = list(g = matrix(.observed_sd_yg^2),
                                             theta = .gp_mod_g$theta)))
      
    } else if (stringr::str_detect(.gp_mod_temp_g$msg, "ERROR")){
      .gp_mod_g <- suppressWarnings(mleHomGP(.rescaled_model_data_g,
                            .rescaled_model_data_g$Z,
                            lower = rep(sqrt(.Machine$double.eps), ncol(.matrix_x)),
                            upper = rep(sqrt(ncol(.matrix_x)), ncol(.matrix_x)),
                            known = list(g = .gp_mod_g$g,
                                         theta = .gp_mod_g$theta)))
                                         #beta0 = .gp_mod_g$beta0))
    } else {
      .gp_mod_g <- .gp_mod_temp_g
    }

    
if(debug == TRUE){
  print("made it past additional model fitting") 
}
   

    #------------- Obtain best x and f ------------------------#
    #------------- subject to satisfying constraint on g ------#
    
    .predictions_f_at_x <- predict(.gp_mod_f,
                                   x = as.matrix(.new_x_grid %>% select(-constraint_threshold_g))) 

    .predictions_g_at_x <- predict(.gp_mod_g,
                                   x = as.matrix(.new_x_grid %>% select(-constraint_threshold_g))) 

    
if(debug == TRUE){
  print("made it up to minimizer")
}

    
    .minimizer <- .new_x_grid %>%
      # convert back to original scale
      mutate(posterior_mean_f_at_best_x = .predictions_f_at_x$mean*.observed_sd_yf + .observed_mean_yf,
             posterior_sd2_f_at_best_x = .predictions_f_at_x$sd2*(.observed_sd_yf)^2, 
             posterior_mean_g_at_best_x = .predictions_g_at_x$mean*.observed_sd_yg + .observed_mean_yg,
             posterior_sd2_g_at_best_x = .predictions_g_at_x$sd2*(.observed_sd_yg)^2) %>%
      mutate(probability_satisfying_constraint_g = pmap_dbl(list(.$posterior_mean_g_at_best_x,
                                                                 sqrt(.$posterior_sd2_g_at_best_x),
                                                                 .$constraint_threshold_g),
                                                            function(x,y,z) pnorm(q = z,
                                                                                  mean = x,
                                                                                  sd = y)),
             feasible_best_x = if_else(probability_satisfying_constraint_g > 
                                    constraint_probability_threshold, 1, 0)) %>%
      group_by(!!!.z_columns_names_syms)

    
if(debug == TRUE){
  print("made it past minimizer")  
}
    
  
   
    .feasible_x_check <- .minimizer %>%
      summarize(no_feasible_x_in_stratum = if_else(all(feasible_best_x == 0), 1L, 0L)) 
    
    .minimizer_temp <- .minimizer %>%
      left_join(.feasible_x_check, by = .feasible_x_check %>% select(-no_feasible_x_in_stratum) %>% names())

    if(all(.feasible_x_check$no_feasible_x_in_stratum == 1L)){
      .minimizer <- .minimizer_temp %>%
        filter(probability_satisfying_constraint_g == max(probability_satisfying_constraint_g)) %>%
      slice_sample(n = 1) 

    } else if(all(.feasible_x_check$no_feasible_x_in_stratum == 0L)){
      .minimizer <- .minimizer_temp %>%
        filter(feasible_best_x == 1) %>%
        filter(posterior_mean_f_at_best_x == min(posterior_mean_f_at_best_x)) %>%
        slice_sample(n = 1)

    } else {
      .minimizer_feasible <-  .minimizer_temp %>%
        filter(no_feasible_x_in_stratum == 0) %>%
        filter(feasible_best_x == 1) %>%
        filter(posterior_mean_f_at_best_x == min(posterior_mean_f_at_best_x)) %>%
        slice_sample(n = 1)
      
      .minimizer_no_feasible <-  .minimizer_temp %>%
        filter(no_feasible_x_in_stratum == 1) %>%
        filter(probability_satisfying_constraint_g == max(probability_satisfying_constraint_g)) %>%
        slice_sample(n = 1)
      
      .minimizer <- .minimizer_feasible %>%
        bind_rows(.minimizer_no_feasible) #%>%
        
    }

    
    
if(debug == TRUE){
  print("made it past three part section")  
}
  
    
    .minimizer <- .minimizer %>%
      ungroup() %>%
      select(-any_of("\"DROP_COLUMN\"")) %>%
      rename_with(~ paste0("best_", .x_columns_names), starts_with("x"))
    .x_columns <- .minimizer %>% select(starts_with("best_x"))
    .minimizer <- .minimizer %>%
      mutate(best_x = .x_columns %>% transpose() %>% map(~ unlist(.), use.names = TRUE))

    if(algorithm == "standard"){

      .strata <- distribution_parameters %>% select(num_range("z",1:20))
      .minimizer <- .minimizer %>%
        bind_cols(.strata)
      
    }

    
    
if(debug == TRUE){
  print("made it up to new_information")
}

    .new_information <- .minimizer %>%
                        right_join(.new_data, by = c(names(.z_columns), "constraint_threshold_g"))
     
    .new_information <- .new_information %>%
      mutate(gp_message_f = .gp_mod_f$msg,
             gp_message_g = .gp_mod_g$msg,
             true_f_at_best_x = map2_dbl(.$best_x, .$z, objective_function, distribution_parameters),
             true_g_at_best_x = map2_dbl(.$best_x, .$z, constraint_function, distribution_parameters))
    
if(debug == TRUE){
  print("made it past new_information")
}

    #--------------- Obtain performance metrics ------------------------#
    .performance <- foreach(i = 1:nrow(.new_information), .combine = 'rbind') %do% {
      quantify_performance(z = .new_information$z[i],
                           best_x = .new_information$best_x[i],
                           mean_estimate_f = .new_information$posterior_mean_f_at_best_x[i],
                           sd2_estimate_f = .new_information$posterior_sd2_f_at_best_x[i],
                           mean_estimate_g = .new_information$posterior_mean_g_at_best_x[i],
                           sd2_estimate_g = .new_information$posterior_sd2_g_at_best_x[i],
                           true_f_at_best_x = .new_information$true_f_at_best_x[i],
                           true_g_at_best_x = .new_information$true_g_at_best_x[i],
                           distribution_parameters = distribution_parameters,
                           x_grid_granularity = x_grid_granularity,
                           seed = seed)
    }

    
    .new_information <- bind_cols(.new_information,
                                 .performance)
    

    .evaluated_data <- bind_rows(.evaluated_data,
                                 .new_information)

    
    
if(debug == TRUE){
  print(paste0("end iteration: ", .iteration)) 
}

if(debug == TRUE){
  if(.iteration == 50)
    return(.evaluated_data)
}   
    
    #--------- Update algorithm stopping criteria -----------------------#
    .n_enrolled <- nrow(.evaluated_data %>% filter(!is.na(f))) #nrow(.evaluated_data %>% filter(af_early_stop == 0, constraint_early_stop == 0))
    .all_strata_stopped <- .evaluated_data %>% 
      group_by(!!!.z_columns_names_syms) %>%
      slice(n()) %>% ungroup() %>% select(af_early_stop, constraint_early_stop) %>% mutate(any_stop = if_else(af_early_stop == 1 | constraint_early_stop == 1, 1, 0)) #no_feasible_x_in_stratum
    .all_strata_stopped_logical <- all(.all_strata_stopped$any_stop == 1) #all(.all_strata_stopped$af_early_stop == 1) | all(.all_strata_stopped$constraint_early_stop == 1) #.all_strata_stopped$no_feasible_x_in_stratum 
    .iteration <- .iteration + 1

if(debug == TRUE){
  print(paste0("n_enrolled: ", .n_enrolled, ", all_strata_stopped_logical: ", .all_strata_stopped_logical)) 
}  
    
    # in case while loop getting stuck
    if(.iteration > max_sample_size * reps_each_x){
      break
    }
    
  }
    
  # Return final dataset
  .evaluated_data %>%
    mutate(simulation = simulation,
           algorithm = algorithm,
           seed = seed,
           true_noise_sd_f = true_noise_sd_f,
           true_noise_sd_g = true_noise_sd_g,
           initial_n_doses = initial_n_doses,
           x_grid_granularity = x_grid_granularity,
           reps_each_x = reps_each_x,
           max_sample_size = max_sample_size,
           acquisition_function = acquisition_function,
           x_generation_type = x_generation_type,
           z_enrollment_type = z_enrollment_type,
           af_stop_threshold = af_stop_threshold,
           constraint_probability_threshold = constraint_probability_threshold,
           rho = rho,
           evaluated_x_violates_constraint = if_else(g > constraint_threshold_g, 1, 0)) 
  
}

