library(RMTIF)

# Known two-sample projection variance: use subject means rather than treating
# all pairwise observations as independent.
net <- matrix(c(1, 2, 3, 4), nrow = 2L, byrow = TRUE)
fit <- rmtif_projection_inference(net, conf_level = 0.95)
stopifnot(isTRUE(all.equal(fit$estimate, 2.5)),
          isTRUE(all.equal(fit$se, sqrt(1.25))),
          length(fit$conf_int) == 2L,
          fit$p_value > 0, fit$p_value < 1)

# Right continuity, bounded time units, and terminal-state continuation.
control <- rmtif_subject(0, c(0, 1.5), c(2, 3), tau = 3.7,
                        death_time = 1.5, death_state = 3, id = "c1")
treatment <- rmtif_subject(1, c(0, 2.1), c(1, 2), tau = 3.7,
                          id = "t1")
stopifnot(identical(unname(rmtif_state_at(control, c(0, 1.5, 3))),
                    c(2, 3, 3)))
pair <- rmtif_pair_times(control, treatment, tau = 3.7, censoring = "none")
stopifnot(isTRUE(all.equal(unname(pair["win"]), 3.7)),
          isTRUE(all.equal(unname(pair["loss"]), 0)))

# Finite-state state-count evaluation and the direct pair estimator agree.
subjects <- list(control, control, treatment, treatment)
subjects[[2L]]$id <- "c2"
subjects[[4L]]$id <- "t2"
direct <- rmtif_complete(subjects, tau = 3.7, conf_level = 0.95)
curve <- rmtif_complete_state_count_curve(subjects, times = c(1, 3.7),
                                         conf_level = 0.95)
stopifnot(isTRUE(all.equal(direct$estimate, 3.7)),
          isTRUE(all.equal(tail(curve$estimate, 1L), direct$estimate)))

# Without censoring the terminal-aware IPCW point estimator agrees with the
# complete-path estimator, and all three computational forms agree.
equivalence <- rmtif_ipcw_equivalence_check(subjects, tau = 3.7,
                                            minimum_survival = 0,
                                            tolerance = 1e-10)
stopifnot(equivalence$absolute_difference < 1e-10,
          equivalence$cached_absolute_difference < 1e-10,
          equivalence$bootstrap_cached_absolute_difference < 1e-10)
ipcw <- rmtif_ipcw(subjects, tau = 3.7, bootstrap = 5L,
                   conf_level = 0.95, bootstrap_success_fraction = 1,
                   seed = 379L)
stopifnot(isTRUE(all.equal(ipcw$estimate, direct$estimate)),
          isTRUE(all.equal(ipcw$se, 0)))

# Exercise refitted arm-specific censoring weights with an observed death
# before censoring and independent losses in both arms.
observed <- list(
  control,
  rmtif_subject(0, c(0, 1), c(2, 1), tau = 3.7,
                censor_time = 2, censor_event = TRUE, id = "c2"),
  rmtif_subject(0, 0, 2, tau = 3.7, id = "c3"),
  treatment,
  rmtif_subject(1, 0, 1, tau = 3.7,
                censor_time = 1.6, censor_event = TRUE, id = "t2"),
  rmtif_subject(1, c(0, 2.5), c(1, 2), tau = 3.7, id = "t3"))
weighted <- rmtif_ipcw_equivalence_check(observed, tau = 3.7,
                                        tolerance = 1e-10)
stopifnot(weighted$absolute_difference < 1e-10,
          weighted$cached_absolute_difference < 1e-10,
          weighted$bootstrap_cached_absolute_difference < 1e-10)

# Both transition specifications produce valid generators. A symmetric
# two-arm generator has zero population net-time contrast.
rates <- matrix(c(0, 0.17, 0.29, 0), nrow = 2L, byrow = TRUE)
generator <- rmtif_ctmc_generator(rates)
stopifnot(isTRUE(all.equal(unname(rowSums(generator)), c(0, 0))),
          isTRUE(all.equal(unname(generator[1L, 2L]), 0.17)))
transition <- matrix(c(0.37, 0.63, 0.44, 0.56), nrow = 2L, byrow = TRUE)
update <- rmtif_poisson_update_generator(transition, rates = c(1.3, 2.1))
stopifnot(isTRUE(all.equal(unname(rowSums(update)), c(0, 0))),
          isTRUE(all.equal(unname(update[1L, 2L]), 1.3 * 0.63)))
truth <- rmtif_true_ctmc_contrast(generator, generator, initial_state = 0,
                                 tau = 3.7, grid_points = 51L)
stopifnot(abs(truth) < 1e-12)

# An explicit simulation seed is reproducible and preserves the caller's RNG.
set.seed(174L)
caller_rng <- .Random.seed
sim1 <- rmtif_simulate_gumbel_trial(n_per_arm = 12L,
  beta_true = c(0.2, 0.3), tau = 3.7, sigma_x = matrix(1), seed = 397L)
stopifnot(identical(.Random.seed, caller_rng))
sim2 <- rmtif_simulate_gumbel_trial(n_per_arm = 12L,
  beta_true = c(0.2, 0.3), tau = 3.7, sigma_x = matrix(1), seed = 397L)
stopifnot(identical(sim1, sim2), length(sim1$subjects) == 24L)

# Zero exit rates yield undefined next-transition probabilities.
# Response confirmation at the final administrative assessment remains valid
# when it is strictly before progression or death.
transition <- data.frame(from = 1, to = 4, baseline_rate = 0,
                          covariate_log_rate = 0, treatment_log_rate = 0)
design <- rmtif_oncology_design(transition,
  states = c(CR = 0L, PR = 1L, SD = 2L, PD = 3L, death = 4L),
  initial_state = 1L, death_state = 4L, scenario = "boundary-test")
intensity <- rmtif_oncology_transition_rates(design, arm = 1, covariate = 0)
stopifnot(intensity$intensity == 0,
          is.na(intensity$next_transition_probability))
trial <- rmtif_simulate_oncology_trial(
  n_per_arm = 2L, design = design, tau = 1.6,
  assessment_interval = 0.8, assessment_jitter = 0, censoring_rate = 0,
  covariate_mean = 0, covariate_sd = 1, seed = 397L)
stopifnot(all(trial$data$objective_response == 1L),
          all(trial$data$confirmed_disease_control == 1L),
          all(trial$data$pfs_event == 0L), all(trial$data$os_event == 0L))

# Public function validation rejects an invalid confidence level.
invalid <- try(rmtif_projection_inference(net, conf_level = 1), silent = TRUE)
stopifnot(inherits(invalid, "try-error"))
