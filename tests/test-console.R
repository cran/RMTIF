library(RMTIF)

# Study functions return structured results without unsolicited console output,
# including checkpoint writes and the startup of parallel scenario workers.
if (requireNamespace("MASS", quietly = TRUE) &&
    requireNamespace("survival", quietly = TRUE)) (function() {
  transition <- data.frame(
    from = c(2, 2, 2, 1, 1, 1, 0, 0, 3),
    to = c(1, 3, 4, 0, 3, 4, 3, 4, 4),
    baseline_rate = c(.23, .17, .07, .21, .16, .08, .13, .06, .27),
    covariate_log_rate = c(-.22, .31, .18, -.22, .31, .18, .31, .18, .18),
    treatment_log_rate = rep(0, 9))
  scenarios <- c("quiet-a", "quiet-b")
  designs <- lapply(scenarios, function(scenario) rmtif_oncology_design(
    transition, states = c(CR = 0L, PR = 1L, SD = 2L, PD = 3L, death = 4L),
    initial_state = 2L, death_state = 4L, scenario = scenario))
  names(designs) <- scenarios
  settings <- list(scenarios = scenarios, designs = designs,
    n_per_arm = 12L, replications = 1L, tau = 3.7,
    assessment_interval = .8, assessment_jitter = .03,
    censoring_rate = 0, covariate_mean = .15, covariate_sd = .9,
    rmt_inference = "complete_projection", regression_bootstrap = 0L,
    seed = 397L, conf_level = .91, bootstrap_success_fraction = .8,
    minimum_survival = 0, test_level = .09)
  checkpoint_dir <- tempfile("rmtif-quiet-checkpoints-")
  on.exit(unlink(checkpoint_dir, recursive = TRUE), add = TRUE)
  serial_text <- capture.output(serial <- do.call(rmtif_oncology_study,
    c(settings, list(workers = 1L, checkpoint_dir = checkpoint_dir))))
  stopifnot(!length(serial_text),
    nrow(serial$summary) == length(scenarios) * 6L,
    nrow(serial$replicates) == length(scenarios) * 6L,
    length(list.files(checkpoint_dir, pattern = "[.]rds$")) == 2L)
  resumed_text <- capture.output(resumed <- do.call(rmtif_oncology_study,
    c(settings, list(workers = 1L, checkpoint_dir = checkpoint_dir))))
  stopifnot(!length(resumed_text),
    identical(serial$summary, resumed$summary),
    identical(serial$replicates, resumed$replicates),
    all(resumed$runtime$resumed_after[1:2] == settings$replications))
  parallel_text <- capture.output(parallel_result <- do.call(
    rmtif_oncology_study, c(settings, list(workers = 2L))))
  fields <- setdiff(names(serial$replicates), "trial_elapsed_seconds")
  stopifnot(!length(parallel_text),
    identical(serial$summary, parallel_result$summary),
    identical(serial$replicates[fields], parallel_result$replicates[fields]))
})()
