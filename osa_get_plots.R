## Set of functions for plotting both regular/personalized algorithms

pacman::p_load(dplyr, tidyr, purrr, ggplot2, cowplot, foreach, doParallel, latex2exp, cli)

##################### Data manipulation #############################

sim_results_constraint_1.5_2_standard <- readRDS("2024-04-05/Data/rwa_early_stopping_toxthresh_1.5_2_standard_rho_0.5.RDS")
sim_results_constraint_1.5_2_personalized <- readRDS("2024-04-05/Data/rwa_early_stopping_toxthresh_1.5_2_personalized_rho_0.5.RDS")

##---------- Extract the results versus the errors -----------#
results_standard <- sim_results_constraint_1.5_2_standard[,1] %>%
  list_rbind() 
results_personalized <- sim_results_constraint_1.5_2_personalized[,1] %>%
  list_rbind() 

#-----------------------------------------------------------------------------#

results_no_missing_standard <- results_standard %>% filter(!is.na(f), algorithm == "standard") %>%
  group_by(algorithm, seed, reps_each_x, af_stop_threshold, true_noise_sd_f, true_noise_sd_g) %>%
  mutate(evaluation = row_number()) %>%
  filter(evaluation == max(evaluation)) %>%
  ungroup() %>%
  mutate(z1 = factor(z1),
         true_noise_sd_f = factor(true_noise_sd_f),
         true_noise_sd_g = factor(true_noise_sd_g),
         constraint_threshold_g = factor(constraint_threshold_g),
         reps_each_x = factor(reps_each_x), 
         af_stop_threshold = factor(af_stop_threshold),
         n_stop = factor(case_when(reps_each_x == 2 & af_stop_threshold == 0.077 ~ 40,
                                   reps_each_x == 2 & af_stop_threshold == 0.060 ~ 60,
                                   reps_each_x == 2 & af_stop_threshold == 0.050 ~ 80,
                                   reps_each_x == 4 & af_stop_threshold == 0.100 ~ 40,
                                   reps_each_x == 4 & af_stop_threshold == 0.073 ~ 60,
                                   reps_each_x == 4 & af_stop_threshold == 0.055 ~ 80)),
         design = factor(if_else(reps_each_x == 2, "S2", "S4")))


results_no_missing_personalized <- results_personalized %>% filter(!is.na(f), algorithm == "personalized") %>%
  group_by(algorithm, z1, seed, reps_each_x, af_stop_threshold, true_noise_sd_f, true_noise_sd_g) %>%
  mutate(evaluation = row_number()) %>%
  filter(evaluation == max(evaluation)) %>%
  ungroup() %>%
  mutate(z1 = factor(z1),
         true_noise_sd_f = factor(true_noise_sd_f),
         true_noise_sd_g = factor(true_noise_sd_g),
         constraint_threshold_g = factor(constraint_threshold_g),
         reps_each_x = factor(reps_each_x), 
         af_stop_threshold = factor(af_stop_threshold),
         n_stop = factor(case_when(reps_each_x == 1 & af_stop_threshold == 0.100 ~ 40,
                                   reps_each_x == 1 & af_stop_threshold == 0.080 ~ 60,
                                   reps_each_x == 1 & af_stop_threshold == 0.065 ~ 80,
                                   reps_each_x == 2 & af_stop_threshold == 0.137 ~ 40,
                                   reps_each_x == 2 & af_stop_threshold == 0.087 ~ 60,
                                   reps_each_x == 2 & af_stop_threshold == 0.070 ~ 80)),
         design = factor(if_else(reps_each_x == 1, "P1", "P2")))


#------------------- expected sample size ----------------------------------------------------------#

ess_standard <- results_no_missing_standard %>%
  group_by(algorithm, reps_each_x, af_stop_threshold, true_noise_sd_f, true_noise_sd_g, 
           constraint_threshold_g, n_stop, design) %>%
  summarize(expected_sample_size = mean(evaluation), .groups = "drop") 

ess_personalized <- results_no_missing_personalized %>%
  group_by(algorithm, z1, reps_each_x, af_stop_threshold, true_noise_sd_f, true_noise_sd_g, 
           constraint_threshold_g, n_stop, design) %>%
  summarize(expected_sample_size = mean(evaluation), .groups = "drop") 

summarized_data_standard <- results_no_missing_standard %>%
  #filter(!is.na(f)) %>% #because this denotes early stopping
  group_by(algorithm, z1, reps_each_x, af_stop_threshold,
           true_noise_sd_f, true_noise_sd_g,  
           constraint_threshold_g, n_stop, design) %>%
  summarise(mean_rpsel_true_f_at_optimal = mean(rpsel_true_f_at_optimal, na.rm = TRUE),
            mean_rpsel_true_g_at_optimal = mean(rpsel_true_g_at_optimal, na.rm = TRUE),
            mean_rpsel_true_f_at_best_x = mean(rpsel_true_f_at_best_x, na.rm = TRUE),
            mean_rpsel_true_g_at_best_x = mean(rpsel_true_g_at_best_x, na.rm = TRUE),
            mean_f = mean(posterior_mean_f_at_best_x, na.rm = TRUE),
            mean_g = mean(posterior_mean_g_at_best_x, na.rm = TRUE),
            mean_dose_units_from_optimal = mean(dose_units_from_optimal, na.rm = TRUE),
            mean_absolute_deviation_f = mean(abs(posterior_mean_f_at_best_x - true_f_at_best_x), na.rm = TRUE),
            mean_absolute_deviation_g = mean(abs(posterior_mean_g_at_best_x - true_g_at_best_x), na.rm = TRUE),
            mean_max_af_value = mean(max_af_value, na.rm = TRUE),
            median_max_af_value = median(max_af_value, na.rm = TRUE),
            .groups = "drop") 
levels(summarized_data_standard$z1) <- c(`0` = TeX("$z_1 = 0$"), `1` = TeX("$z_1 = 1$"))





summarized_data_personalized <- results_no_missing_personalized %>%
  group_by(algorithm, z1, reps_each_x, af_stop_threshold,
           true_noise_sd_f, true_noise_sd_g,  
           constraint_threshold_g, n_stop, design) %>%
  summarise(mean_rpsel_true_f_at_optimal = mean(rpsel_true_f_at_optimal, na.rm = TRUE),
            mean_rpsel_true_g_at_optimal = mean(rpsel_true_g_at_optimal, na.rm = TRUE),
            mean_rpsel_true_f_at_best_x = mean(rpsel_true_f_at_best_x, na.rm = TRUE),
            mean_rpsel_true_g_at_best_x = mean(rpsel_true_g_at_best_x, na.rm = TRUE),
            mean_f = mean(posterior_mean_f_at_best_x, na.rm = TRUE),
            mean_g = mean(posterior_mean_g_at_best_x, na.rm = TRUE),
            mean_dose_units_from_optimal = mean(dose_units_from_optimal, na.rm = TRUE),
            mean_absolute_deviation_f = mean(abs(posterior_mean_f_at_best_x - true_f_at_best_x), na.rm = TRUE),
            mean_absolute_deviation_g = mean(abs(posterior_mean_g_at_best_x - true_g_at_best_x), na.rm = TRUE),
            mean_max_af_value = mean(max_af_value, na.rm = TRUE),
            median_max_af_value = median(max_af_value, na.rm = TRUE),
            .groups = "drop")
levels(summarized_data_personalized$z1) <- c(`0` = TeX("$z_1 = 0$"), `1` = TeX("$z_1 = 1$"))


summarized_data <- summarized_data_standard %>%
  bind_rows(summarized_data_personalized) %>%
  mutate(design = factor(design, levels = c("P1", "P2", "S2", "S4")))


#-------------- Plots ---------------------------------------------------------#
#----------- Dose units from optimal ------------------------------------------#

plot_exp_dose_units <- ggplot(summarized_data, 
       aes(x = n_stop, y = mean_dose_units_from_optimal, 
           group = design,
           color = design,
           linetype = design,
           shape = constraint_threshold_g)) + 
  facet_grid(z1 ~ ., labeller = label_parsed, scales = "free") +
  geom_point(size = 1.7) + 
  geom_line() +
  scale_y_continuous(limits = c(NA,NA)) +
  scale_color_manual("algorithm", values = c("#26495c", "#468B97", "#c66b3d", "#EF6262")) +
  scale_linetype_manual(name = "algorithm",
                        values = c("F1", "dotdash", "dashed", "dotted")) +
  scale_shape_manual(parse(text = TeX("$g^{\\dagger}$")),
                     values = c(0,6),
                     labels = c("1.5", "2")) +
  theme_bw() +
  guides(color = guide_legend(nrow = 4, order = 1, override.aes=list(linewidth = 1, shape = NA)),
         linetype = guide_legend(nrow = 4, order = 1)) +
  ylab("dose units from optimal") +
  xlab(parse(text = TeX("$n_{STOP}$"))) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        #legend.title = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 7),
        aspect.ratio = 0.67,
        plot.margin = unit(c(0.5, 0.25, 0, 0.25), "cm"))


#------------------ RPSEL at optimal ---------------------------------#

plot_rpsel <- ggplot(summarized_data, 
       aes(x = n_stop, y = mean_rpsel_true_f_at_optimal, 
           group = design,
           color = design,
           linetype = design,
           shape = constraint_threshold_g)) + 
  facet_grid(z1 ~ ., labeller = label_parsed, scales = "free") +
  geom_point(size = 1.7) + 
  geom_line() +
  scale_y_continuous(limits = c(NA,NA)) +
  scale_color_manual("algorithm", values = c("#26495c", "#468B97", "#c66b3d", "#EF6262")) +
  scale_linetype_manual(name = "algorithm",
                        values = c("F1", "dotdash", "dashed", "dotted")) +
  scale_shape_manual(parse(text = TeX("$g^{\\dagger}$")),
                     values = c(0,6),
                     labels = c("1.5", "2")) +
  theme_bw() +
  guides(color = guide_legend(nrow = 1, order = 1, override.aes=list(linewidth = 1, shape = NA)),
         linetype = guide_legend(nrow = 1, order = 1)) +
  ylab("rpsel") +
  xlab(parse(text = TeX("$n_{STOP}$"))) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        legend.position = "none",
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 7),
        aspect.ratio = 0.67,
        plot.margin = unit(c(0.5, 0.25, 0, 0.25), "cm"))


### get legend for plot below
### change legend.position above to "top"
set_null_device(cairo_ps) 
plot_legend <- get_plot_component(ggplot(summarized_data, 
                                         aes(x = n_stop, y = mean_rpsel_true_f_at_optimal, 
                                             group = design,
                                             color = design,
                                             linetype = design,
                                             shape = constraint_threshold_g)) + 
                                    facet_grid(z1 ~ ., labeller = label_parsed, scales = "free") +
                                    geom_point(size = 1.7) + 
                                    geom_line() +
                                    scale_y_continuous(limits = c(NA,NA)) +
                                    scale_color_manual("algorithm", values = c("#26495c", "#468B97", "#c66b3d", "#EF6262")) +
                                    scale_linetype_manual(name = "algorithm",
                                                          values = c("F1", "dotdash", "dashed", "dotted")) +
                                    scale_shape_manual(parse(text = TeX("$g^{\\dagger}$")),
                                                       values = c(0,6),
                                                       labels = c("1.5", "2")) +
                                    theme_bw() +
                                    guides(color = guide_legend(nrow = 1, order = 1, override.aes=list(linewidth = 1, shape = NA)),
                                           linetype = guide_legend(nrow = 1, order = 1)) +
                                    ylab("rpsel") +
                                    xlab(parse(text = TeX("$n_{STOP}$"))) +
                                    theme(panel.grid.minor = element_blank(),
                                          panel.grid.major.x = element_blank(),
                                          legend.position = "top",
                                          axis.title = element_text(size = 10),
                                          axis.text = element_text(size = 7),
                                          aspect.ratio = 0.67,
                                          plot.margin = unit(c(0,0,0,0), "cm")),
                                  'guide-box', return_all = TRUE)[[4]]

ggdraw(plot_legend)

#-------------- Expected Sample Size-------------------------------------------#

ess_all <- ess_standard %>%
  bind_rows(ess_personalized) %>%
  mutate(design = factor(design, levels = c("P1", "P2", "S2", "S4")))
levels(ess_all$n_stop) <- c(`40` = TeX("$n_{STOP}=40$"),
                                    `60` = TeX("$n_{STOP}=60$"),
                                    `80` = TeX("$n_{STOP}=80$"))

plot_ess <- ggplot(ess_all,
       aes(x = design,
           y = expected_sample_size, 
           fill = z1)) +
  facet_grid(. ~ n_stop, labeller = label_parsed) +
  geom_col() +
  scale_fill_manual(values = c("#79AEB2", "#4A6274"),
                    na.value = "#657580",
                    limits = c("0", "1"),
                    labels = unname(TeX(c("$z_1=0$", "$z_1=1$")))) +
  theme_bw() +
  ylab("expected sample size") +
  xlab("") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        #legend.position = "top",
        legend.title = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        aspect.ratio = 1)

#----------------- Toxic doses ------------------------------------------------#

#------------ Toxic Doses ----------------------------------------------------#

toxic_doses_standard <- results_standard %>% filter(!is.na(f)) %>%
  group_by(algorithm, reps_each_x, af_stop_threshold, seed, true_noise_sd_f, true_noise_sd_g, constraint_threshold_g) %>%
  summarize(sum_evaluated_x_violates_constraint = sum(evaluated_x_violates_constraint, na.rm = TRUE),
            .groups = "keep") %>%
  ungroup(seed) %>%
  summarize(mean_number_parts_toxicity = mean(sum_evaluated_x_violates_constraint),
            .groups = "drop") %>%
  mutate(z1 = NA,
         true_noise_sd_f = factor(true_noise_sd_f),
         true_noise_sd_g = factor(true_noise_sd_g),
         constraint_threshold_g = factor(constraint_threshold_g),
         reps_each_x = factor(reps_each_x), 
         af_stop_threshold = factor(af_stop_threshold),
         n_stop = factor(case_when(reps_each_x == 2 & af_stop_threshold == 0.077 ~ 40,
                                   reps_each_x == 2 & af_stop_threshold == 0.060 ~ 60,
                                   reps_each_x == 2 & af_stop_threshold == 0.050 ~ 80,
                                   reps_each_x == 4 & af_stop_threshold == 0.100 ~ 40,
                                   reps_each_x == 4 & af_stop_threshold == 0.073 ~ 60,
                                   reps_each_x == 4 & af_stop_threshold == 0.055 ~ 80)),
         design = factor(if_else(reps_each_x == 2, "S2", "S4")))

toxic_doses_personalized <- results_personalized %>% filter(!is.na(f)) %>%
  group_by(algorithm, z1, reps_each_x, af_stop_threshold, seed, true_noise_sd_f, true_noise_sd_g, constraint_threshold_g) %>%
  summarize(sum_evaluated_x_violates_constraint = sum(evaluated_x_violates_constraint, na.rm = TRUE),
            .groups = "keep") %>%
  ungroup(seed) %>%
  summarize(mean_number_parts_toxicity = mean(sum_evaluated_x_violates_constraint),
            .groups = "drop") %>%
  mutate(z1 = factor(z1),
         true_noise_sd_f = factor(true_noise_sd_f),
         true_noise_sd_g = factor(true_noise_sd_g),
         constraint_threshold_g = factor(constraint_threshold_g),
         reps_each_x = factor(reps_each_x), 
         af_stop_threshold = factor(af_stop_threshold),
         n_stop = factor(case_when(reps_each_x == 1 & af_stop_threshold == 0.100 ~ 40,
                                   reps_each_x == 1 & af_stop_threshold == 0.080 ~ 60,
                                   reps_each_x == 1 & af_stop_threshold == 0.065 ~ 80,
                                   reps_each_x == 2 & af_stop_threshold == 0.137 ~ 40,
                                   reps_each_x == 2 & af_stop_threshold == 0.087 ~ 60,
                                   reps_each_x == 2 & af_stop_threshold == 0.070 ~ 80)),
         design = factor(if_else(reps_each_x == 1, "P1", "P2")))

toxic_doses_all <- toxic_doses_personalized %>%
  bind_rows(toxic_doses_standard)
levels(toxic_doses_all$n_stop) <- c(`40` = TeX("$n_{STOP}=40$"),
                                    `60` = TeX("$n_{STOP}=60$"),
                                    `80` = TeX("$n_{STOP}=80$"))

plot_toxic_doses <- ggplot(toxic_doses_all,
                           aes(x = design,
                               y = mean_number_parts_toxicity,
  fill = factor(z1))) +
  facet_grid(. ~ n_stop, labeller = label_parsed) +
  geom_col() +
  scale_fill_manual(values = c("#79AEB2", "#4A6274"),
                    na.value = "#657580",
                    limits = c("0", "1"),
                    labels = unname(TeX(c("$z_1=0$", "$z_1=1$")))) +
  theme_bw() +
  ylab("toxic doses") +
  xlab("") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        #legend.position = "top",
        legend.title = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        aspect.ratio = 1)


#-------------- Prob incorrectly determining no feasible doses ----------------#
#-------------- AKA stopping early for toxicity ----------------#


stopped_tox_standard <- results_standard %>%
  filter(constraint_early_stop == 1) %>%
  select(algorithm, reps_each_x, seed, simulation, af_stop_threshold, true_noise_sd_f, true_noise_sd_g) %>%
  distinct() %>%
  group_by(algorithm, reps_each_x, af_stop_threshold, true_noise_sd_f, true_noise_sd_g) %>%
  select(group_cols()) %>%
  tally() %>%
  mutate(prop = n/1100,
         true_noise_sd_f = factor(true_noise_sd_f),
         true_noise_sd_g = factor(true_noise_sd_g),
         reps_each_x = factor(reps_each_x), 
         af_stop_threshold = factor(af_stop_threshold),
         n_stop = factor(case_when(reps_each_x == 2 & af_stop_threshold == 0.077 ~ 40,
                                   reps_each_x == 2 & af_stop_threshold == 0.060 ~ 60,
                                   reps_each_x == 2 & af_stop_threshold == 0.050 ~ 80,
                                   reps_each_x == 4 & af_stop_threshold == 0.100 ~ 40,
                                   reps_each_x == 4 & af_stop_threshold == 0.073 ~ 60,
                                   reps_each_x == 4 & af_stop_threshold == 0.055 ~ 80)),
         design = factor(if_else(reps_each_x == 2, "S2", "S4"))) %>%
  ungroup() %>%
  select(-true_noise_sd_f, -true_noise_sd_g) %>%
  uncount(2) %>%
  mutate(z1 = factor(rep(c(0,1), 6)))

stopped_tox_personalized <- results_personalized %>%
  filter(constraint_early_stop == 1) %>%
  select(algorithm, z1, reps_each_x, seed, simulation, af_stop_threshold, true_noise_sd_f, true_noise_sd_g) %>%
  distinct() %>%
  group_by(algorithm, z1, reps_each_x, af_stop_threshold, true_noise_sd_f, true_noise_sd_g) %>%
  select(group_cols()) %>%
  tally() %>%
  mutate(prop = n/1100,
         z1 = factor(z1),
         reps_each_x = factor(reps_each_x), 
         af_stop_threshold = factor(af_stop_threshold),
         n_stop = factor(case_when(reps_each_x == 1 & af_stop_threshold == 0.100 ~ 40,
                                   reps_each_x == 1 & af_stop_threshold == 0.080 ~ 60,
                                   reps_each_x == 1 & af_stop_threshold == 0.065 ~ 80,
                                   reps_each_x == 2 & af_stop_threshold == 0.137 ~ 40,
                                   reps_each_x == 2 & af_stop_threshold == 0.087 ~ 60,
                                   reps_each_x == 2 & af_stop_threshold == 0.070 ~ 80)),
         design = factor(if_else(reps_each_x == 1, "P1", "P2"))) %>%
  ungroup() %>%
  select(-true_noise_sd_f, -true_noise_sd_g)


stopped_tox <- stopped_tox_personalized %>%
  bind_rows(stopped_tox_standard)
levels(stopped_tox$n_stop) <- c(`40` = TeX("$n_{STOP}=40$"),
                                    `60` = TeX("$n_{STOP}=60$"),
                                    `80` = TeX("$n_{STOP}=80$"))
levels(stopped_tox$z1) <- c(`0` = TeX("$z_1 = 0$"), `1` = TeX("$z_1 = 1$"))
plot_stopped_tox <- ggplot(stopped_tox,
       aes(x = design,
           y = prop,
           fill = design)) +
  facet_grid(z1 ~ n_stop, labeller = label_parsed) +
  geom_col() +
  scale_fill_manual(values = c("#26495c", "#468B97", "#c66b3d", "#EF6262")) +
  theme_bw() +
  ylab("probability no feasible dose") +
  xlab("") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        #legend.position = "none",
        legend.title = element_blank(),
        axis.title = element_text(size = 10),
        axis.text = element_text(size = 10),
        aspect.ratio = 1)



#----------- Unique doses ----------------------------------------------#

results %>% 
  filter(algorithm_name == .algorithm_name) %>%
  select(simulation, algorithm_name, af_stop_threshold, iteration, x1, x2, af_early_stop) %>%
  distinct() %>%
  filter(af_early_stop == 0) %>% #!is.na(x_1_rec)
  select(simulation, algorithm_name, af_stop_threshold, x1, x2) %>%
  distinct() %>%
  select(-x1, -x2) %>%
  group_by(simulation, algorithm_name, af_stop_threshold) %>%
  tally() %>%
  group_by(simulation, algorithm_name, af_stop_threshold) %>%
  summarize(total_doses = sum(n), .groups = "drop") %>% 
  group_by(algorithm_name, af_stop_threshold) %>%
  summarize(mean_doses = mean(total_doses), .groups = "drop") %>%
  mutate(nstop = c(.middle_stop, .early_stop))


avg_number_doses_standard <- results_standard %>% filter(!is.na(f)) %>%
  select(algorithm, reps_each_x, af_stop_threshold, seed, true_noise_sd_f, true_noise_sd_g, x1, x2) %>%
  distinct() %>%
  select(-x1, -x2) %>%
  group_by(algorithm, reps_each_x, af_stop_threshold, seed, true_noise_sd_f, true_noise_sd_g) %>%
  tally() %>%
  summarize(total_doses = sum(n), .groups = "keep") %>%
  ungroup(seed) %>%
  summarize(mean_number_doses = mean(total_doses),
            .groups = "drop") 

# don't consider strata here since same engineering for same doses 
avg_number_doses_personalized <- results_personalized %>% filter(!is.na(f)) %>%
  select(algorithm, reps_each_x, af_stop_threshold, seed, true_noise_sd_f, true_noise_sd_g, x1, x2) %>%
  distinct() %>%
  select(-x1, -x2) %>%
  group_by(algorithm, reps_each_x, af_stop_threshold, seed, true_noise_sd_f, true_noise_sd_g) %>%
  tally() %>%
  summarize(total_doses = sum(n), .groups = "keep") %>%
  ungroup(seed) %>%
  summarize(mean_number_doses = mean(total_doses),
            .groups = "drop") 

toxic_doses_personalized <- results_personalized %>% filter(!is.na(f)) %>%
  group_by(algorithm, z1, reps_each_x, af_stop_threshold, seed, true_noise_sd_f, true_noise_sd_g, constraint_threshold_g) %>%
  summarize(sum_evaluated_x_violates_constraint = sum(evaluated_x_violates_constraint, na.rm = TRUE),
            .groups = "keep") %>%
  ungroup(seed) %>%
  summarize(mean_number_parts_toxicity = mean(sum_evaluated_x_violates_constraint),
            .groups = "drop") %>%
  mutate(z1 = factor(z1),
         true_noise_sd_f = factor(true_noise_sd_f),
         true_noise_sd_g = factor(true_noise_sd_g),
         constraint_threshold_g = factor(constraint_threshold_g),
         reps_each_x = factor(reps_each_x), 
         af_stop_threshold = factor(af_stop_threshold),
         n_stop = factor(case_when(reps_each_x == 1 & af_stop_threshold == 0.100 ~ 40,
                                   reps_each_x == 1 & af_stop_threshold == 0.080 ~ 60,
                                   reps_each_x == 1 & af_stop_threshold == 0.065 ~ 80,
                                   reps_each_x == 2 & af_stop_threshold == 0.137 ~ 40,
                                   reps_each_x == 2 & af_stop_threshold == 0.087 ~ 60,
                                   reps_each_x == 2 & af_stop_threshold == 0.070 ~ 80)),
         design = factor(if_else(reps_each_x == 1, "P1", "P2")))

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
                     x1 = c(0.25, 0.5),
                     x2 = c(0.75, 0.75))
levels(point_data$z1) <- c(`0` = TeX("$z_1 = 0$"), `1` = TeX("$z_1 = 1$"))


## This is plot of data generating function
true_dgm_plot <- ggplot(plot_data, aes(x = x1, y = x2, z = f)) +
  facet_wrap(z1~., labeller = label_parsed) +
  geom_contour_filled(bins = 6) + 
  geomtextpath::geom_textcontour(aes(x = x1, y = x2, z = g),
                                 breaks = c(0.5, 1, 1.5, 2, 2.5, 3),
                                 size = 2.5,
                                 color = "white", 
                                 linetype = "dashed",
                                 inherit.aes = FALSE) +
  geom_point(data = point_data,
             aes(x = x1, y = x2), shape = 8, color = "white", inherit.aes = FALSE) +
  coord_cartesian(xlim = c(-0.1, 1.1), ylim = c(-0.1, 1.1)) +
  guides(fill = guide_legend(title = "", ncol = 1, byrow = FALSE)) +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 2.5, 5, 7.5, 10)) + #match scale of other plots
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = c(0, 25, 50, 75, 100)) +
  scale_fill_viridis_d(option = "D") +
  ylab(parse(text = TeX('$d_2$')))+
  xlab(parse(text = TeX('$d_1$')))+
  theme(panel.background = element_rect(fill='darkgrey'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.text = element_text(size = 8),
        aspect.ratio = 1)

true_dgm_plot

right_col <- plot_grid(plot_ess,
                     plot_toxic_doses,
                      plot_stopped_tox,
                     ncol = 1,
                     rel_heights = c(1,1,1.65),
                     align = "v",
                     axis = "l",
                     labels = c("B", "C", "F"))

left_col <- plot_grid(true_dgm_plot,
                      
                      plot_grid(plot_exp_dose_units,
                                plot_rpsel, nrow = 1, 
                                labels = c("D","E")),
                      plot_legend,
                      rel_heights = c(1,1,0.1),
                      labels = c("A"),
                      ncol = 1)


#---------------------------------------------------------------------------#
set_null_device(cairo_ps) # this is required for the latex
final_plot <- plot_grid(left_col,
                        right_col,
                        ncol = 2,
                        scale = 0.95)

ggsave2(filename = "rwa.eps",
        device = cairo_ps,
        plot = final_plot,
        width = 11,
        height = 8.5,
        units = "in",
        dpi = 800)
