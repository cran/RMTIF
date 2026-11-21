# RMTIF: restricted mean time in favor of treatment.
# Lower numeric state values are always interpreted as better.

# -----------------------------------------------------------------------------
# Validation and subject representation
# -----------------------------------------------------------------------------

.rmtif_scalar_numeric <- function(x, name, lower = -Inf, finite = TRUE) {
  if (!is.numeric(x) || length(x) != 1L || is.na(x) ||
      (finite && !is.finite(x)) || x < lower) {
    stop(name, " must be a numeric scalar >= ", lower, ".", call. = FALSE)
  }
  as.numeric(x)
}

.rmtif_scalar_flag <- function(x, name) {
  if (!is.logical(x) || length(x) != 1L || is.na(x)) {
    stop(name, " must be TRUE or FALSE.", call. = FALSE)
  }
  x
}

.rmtif_scalar_count <- function(x, name, lower = 0L) {
  x <- .rmtif_scalar_numeric(x, name, lower = lower)
  if (x != floor(x) || x > .Machine$integer.max) {
    stop(name, " must be a representable whole number.", call. = FALSE)
  }
  as.integer(x)
}

.rmtif_conf_level <- function(x) {
  x <- .rmtif_scalar_numeric(x, "conf_level", lower = 0)
  if (x <= 0 || x >= 1) {
    stop("conf_level must be strictly between 0 and 1.", call. = FALSE)
  }
  x
}

.rmtif_bootstrap_success_fraction <- function(x) {
  x <- .rmtif_scalar_numeric(x, "bootstrap_success_fraction", lower = 0)
  if (x <= 0 || x > 1) {
    stop("bootstrap_success_fraction must be in (0, 1].", call. = FALSE)
  }
  x
}

.rmtif_minimum_survival <- function(x) {
  x <- .rmtif_scalar_numeric(x, "minimum_survival", lower = 0)
  if (x > 1) {
    stop("minimum_survival must be between 0 and 1.", call. = FALSE)
  }
  x
}

.rmtif_condition_number <- function(x) {
  singular_values <- svd(as.matrix(x), nu = 0L, nv = 0L)$d
  if (!length(singular_values) || any(!is.finite(singular_values))) {
    return(Inf)
  }
  if (min(singular_values) <= sqrt(.Machine$double.eps) * max(singular_values)) {
    return(Inf)
  }
  max(singular_values) / min(singular_values)
}

.rmtif_pinv <- function(x, tolerance = sqrt(.Machine$double.eps)) {
  x <- as.matrix(x)
  z <- svd(x)
  if (!length(z$d)) {
    stop("Cannot invert an empty matrix.", call. = FALSE)
  }
  cutoff <- tolerance * max(z$d, 1)
  inverse <- ifelse(z$d > cutoff, 1 / z$d, 0)
  z$v %*% (inverse * t(z$u))
}

.rmtif_rng_spec <- function() {
  c(
    kind = "Mersenne-Twister",
    normal.kind = "Inversion",
    sample.kind = "Rejection"
  )
}

.rmtif_set_rng_spec <- function(spec = .rmtif_rng_spec()) {
  current <- base::RNGkind()
  if (length(current) < 3L) {
    stop(
      "The reproducible simulation streams require R >= 3.6.0 so that ",
      "sample.kind can be fixed explicitly.", call. = FALSE
    )
  }
  required <- c("kind", "normal.kind", "sample.kind")
  if (!is.character(spec) || !all(required %in% names(spec)) ||
      any(!nzchar(spec[required]))) {
    stop("spec must name kind, normal.kind, and sample.kind.", call. = FALSE)
  }
  do.call(base::RNGkind, as.list(unname(spec[required])))
  invisible(spec[required])
}

.rmtif_with_seed <- function(seed, code) {
  if (is.null(seed)) {
    return(force(code))
  }
  seed <- .rmtif_scalar_count(seed, "seed", lower = 0L)
  old_kind <- base::RNGkind()
  old_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (old_exists) {
    old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  }
  on.exit({
    do.call(base::RNGkind, as.list(unname(old_kind)))
    if (old_exists) {
      assign(".Random.seed", old_seed, envir = .GlobalEnv)
    } else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
      rm(".Random.seed", envir = .GlobalEnv)
    }
  }, add = TRUE)
  .rmtif_set_rng_spec()
  base::set.seed(seed)
  force(code)
}

.rmtif_next_seed <- function(seed) {
  as.integer((as.double(seed) %% (.Machine$integer.max - 1)) + 1)
}

#' Construct one observed ordered-state trajectory
#'
#' The function creates the common subject representation used by every rmt-IF
#' estimator in this file.  `times` and `states` define a right-continuous,
#' piecewise-constant observed path.  An observed terminal event is inserted at
#' `death_time` and carried forward; it never truncates the outcome path.
#' `censor_time` is loss to follow-up or administrative end, and
#' `censor_event` distinguishes loss to follow-up from administrative end.
#'
#' @param arm Treatment indicator, 0 for control or 1 for treatment.
#' @param times Numeric state-observation/change times and must include baseline.
#' @param states Ordered states corresponding to `times`; smaller is better.
#' @param tau Planned restriction horizon.
#' @param covariates Optional numeric baseline covariate vector.
#' @param censor_time Loss-to-follow-up or administrative time.
#' @param censor_event Whether `censor_time` is a loss-to-follow-up event.
#' @param death_time Observed terminal-event time, or `Inf` if unobserved/absent.
#' @param death_state Ordered absorbing terminal state, or `NULL` if not used.
#' @param id Optional subject identifier.
#' @return An object of class `rmtif_subject`.
#' @export
rmtif_subject <- function(arm, times, states, tau, covariates = numeric(),
                          censor_time = tau, censor_event = FALSE,
                          death_time = Inf, death_state = NULL, id = NULL) {
  tau <- .rmtif_scalar_numeric(tau, "tau", lower = 0)
  if (!(arm %in% c(0, 1)) || length(arm) != 1L || is.na(arm)) {
    stop("arm must be exactly 0 or 1.", call. = FALSE)
  }
  times <- as.numeric(times)
  states <- as.numeric(states)
  if (!length(times) || length(times) != length(states) ||
      anyNA(times) || anyNA(states) || any(!is.finite(times)) ||
      any(!is.finite(states))) {
    stop("times and states must be finite, equally sized, non-empty vectors.",
         call. = FALSE)
  }
  if (any(times < 0) || !any(times == 0)) {
    stop("times must be nonnegative and include baseline time 0.", call. = FALSE)
  }
  censor_time <- .rmtif_scalar_numeric(censor_time, "censor_time", lower = 0,
                                       finite = FALSE)
  censor_event <- .rmtif_scalar_flag(censor_event, "censor_event")
  death_time <- .rmtif_scalar_numeric(death_time, "death_time", lower = 0,
                                      finite = FALSE)
  if (censor_event && censor_time >= tau) {
    censor_event <- FALSE
  }
  if (is.finite(death_time) && death_time > min(censor_time, tau)) {
    stop("death_time must be observed no later than censor_time and tau.",
         call. = FALSE)
  }
  # Once an absorbing terminal event has been observed, a later recorded
  # censoring time is not an observed censoring event. Canonicalizing this
  # flag keeps metadata, censoring KM, and reported LTFU counts consistent.
  if (is.finite(death_time)) {
    censor_event <- FALSE
  }
  if (!censor_event && censor_time < tau && !is.finite(death_time)) {
    stop(
      "censor_time before tau must be marked as a loss-to-follow-up event ",
      "unless an observed death has already ended the path.",
      call. = FALSE
    )
  }
  if (is.finite(death_time) &&
      (is.null(death_state) || length(death_state) != 1L ||
       !is.numeric(death_state) || !is.finite(death_state))) {
    stop("A finite death_time requires one finite numeric death_state.",
         call. = FALSE)
  }
  if (!is.numeric(covariates) || anyNA(covariates) ||
      any(!is.finite(covariates))) {
    stop("covariates must be a finite numeric vector.", call. = FALSE)
  }

  observable_limit <- if (is.finite(death_time)) death_time else min(censor_time, tau)
  keep <- times <= observable_limit & times <= tau
  times <- times[keep]
  states <- states[keep]
  ord <- order(times, seq_along(times))
  times <- times[ord]
  states <- states[ord]
  keep_last <- !duplicated(times, fromLast = TRUE)
  times <- times[keep_last]
  states <- states[keep_last]

  if (is.finite(death_time)) {
    before <- times < death_time
    times <- c(times[before], death_time)
    states <- c(states[before], as.numeric(death_state))
  }
  if (!length(times) || times[1L] != 0) {
    stop("The observed path must contain a state at time 0.", call. = FALSE)
  }

  result <- list(
    id = id,
    arm = as.integer(arm),
    covariates = as.numeric(covariates),
    times = times,
    states = states,
    tau = tau,
    censor_time = censor_time,
    censor_event = censor_event,
    death_time = death_time,
    death_state = if (is.null(death_state)) NA_real_ else as.numeric(death_state)
  )
  class(result) <- "rmtif_subject"
  result
}

.rmtif_validate_subjects <- function(subjects, tau, minimum_per_arm = 1L) {
  tau <- .rmtif_scalar_numeric(tau, "tau", lower = 0)
  if (!is.list(subjects) || !length(subjects) ||
      !all(vapply(subjects, inherits, logical(1), what = "rmtif_subject"))) {
    stop("subjects must be a non-empty list made by rmtif_subject().",
         call. = FALSE)
  }
  arms <- vapply(subjects, `[[`, integer(1), "arm")
  horizons <- vapply(subjects, `[[`, numeric(1), "tau")
  if (any(tau > horizons)) {
    stop("tau exceeds the recorded follow-up horizon for at least one subject.",
         call. = FALSE)
  }
  counts <- tabulate(arms + 1L, nbins = 2L)
  if (any(counts < minimum_per_arm)) {
    stop("At least ", minimum_per_arm, " subject(s) per arm are required.",
         call. = FALSE)
  }
  list(
    control = subjects[arms == 0L],
    treatment = subjects[arms == 1L],
    n0 = counts[1L],
    n1 = counts[2L]
  )
}

#' Evaluate a subject's right-continuous state
#'
#' This helper is intended for transparent checking of the path construction.
#' It returns `NA` after unobserved loss to follow-up, but continues to return
#' the absorbing terminal state after an observed death.
#'
#' @param subject An `rmtif_subject`.
#' @param time One or more nonnegative evaluation times.
#' @return A numeric vector of states, with `NA` where the path is unavailable.
#' @export
rmtif_state_at <- function(subject, time) {
  if (!inherits(subject, "rmtif_subject")) {
    stop("subject must be made by rmtif_subject().", call. = FALSE)
  }
  time <- as.numeric(time)
  if (anyNA(time) || any(time < 0) || any(!is.finite(time))) {
    stop("time must contain finite nonnegative values.", call. = FALSE)
  }
  # findInterval accepts a vector: evaluating a whole observation grid at once
  # avoids one R call per subject per interval without changing the path rule.
  result <- rep(NA_real_, length(time))
  available <- subject$censor_time >= pmin(time, subject$death_time)
  after_death <- available & is.finite(subject$death_time) &
    time >= subject$death_time
  result[after_death] <- subject$death_state
  observed <- available & !after_death
  if (any(observed)) {
    index <- findInterval(time[observed], subject$times)
    result[observed] <- subject$states[index]
  }
  result
}

.rmtif_available <- function(subject, time) {
  subject$censor_time >= pmin(time, subject$death_time)
}

# -----------------------------------------------------------------------------
# Arm-specific censoring Kaplan-Meier and pairwise integration
# -----------------------------------------------------------------------------

.rmtif_censor_km <- function(subjects, tau) {
  followup <- vapply(subjects, function(x) {
    min(x$censor_time, x$death_time, tau)
  }, numeric(1))
  event <- vapply(subjects, function(x) {
    as.integer(x$censor_event &&
                 x$censor_time < min(x$death_time, tau))
  }, integer(1))
  event_times <- sort(unique(followup[event == 1L]))
  if (!length(event_times)) {
    return(structure(list(times = numeric(), survival = numeric(),
                          followup = followup, event = event),
                     class = "rmtif_censor_km"))
  }
  sorted_followup <- sort(followup)
  # With left.open = TRUE, findInterval counts observations strictly below
  # each event time, so the complement is the predictable at-risk count.
  at_risk <- length(followup) - findInterval(
    event_times, sorted_followup, left.open = TRUE
  )
  events <- tabulate(
    match(followup[event == 1L], event_times),
    nbins = length(event_times)
  )
  if (any(at_risk <= 0L)) {
    stop("Empty censoring risk set encountered.", call. = FALSE)
  }
  survival <- cumprod(1 - events / at_risk)
  structure(list(times = event_times, survival = survival,
                 followup = followup, event = event),
            class = "rmtif_censor_km")
}

.rmtif_km_left <- function(km, time) {
  time_names <- names(time)
  time <- as.numeric(time)
  index <- findInterval(time, km$times, left.open = TRUE)
  result <- rep(1, length(index))
  after_first_jump <- index > 0L
  result[after_first_jump] <- km$survival[index[after_first_jump]]
  names(result) <- time_names
  result
}

.rmtif_pair_integral <- function(control, treatment, tau,
                                 km0 = NULL, km1 = NULL,
                                 extra_times = numeric(),
                                 minimum_survival = 0) {
  if (!inherits(control, "rmtif_subject") || control$arm != 0L ||
      !inherits(treatment, "rmtif_subject") || treatment$arm != 1L) {
    stop("The pair must contain an arm-0 control and arm-1 treatment subject.",
         call. = FALSE)
  }
  weighted <- !is.null(km0) || !is.null(km1)
  if (xor(is.null(km0), is.null(km1))) {
    stop("Both arm-specific censoring curves must be supplied together.",
         call. = FALSE)
  }
  grid <- sort(unique(c(
    0, tau,
    control$times[control$times > 0 & control$times < tau],
    treatment$times[treatment$times > 0 & treatment$times < tau],
    control$censor_time[control$censor_event & control$censor_time > 0 &
                          control$censor_time < tau],
    treatment$censor_time[treatment$censor_event &
                            treatment$censor_time > 0 &
                            treatment$censor_time < tau],
    control$death_time[control$death_time > 0 & control$death_time < tau],
    treatment$death_time[treatment$death_time > 0 &
                           treatment$death_time < tau],
    extra_times[extra_times > 0 & extra_times < tau]
  )))
  win <- loss <- 0
  if (length(grid) < 2L) {
    return(c(win = win, loss = loss, net = win - loss))
  }
  for (k in seq_len(length(grid) - 1L)) {
    midpoint <- (grid[k] + grid[k + 1L]) / 2
    width <- grid[k + 1L] - grid[k]
    if (width <= 0 || !.rmtif_available(control, midpoint) ||
        !.rmtif_available(treatment, midpoint)) {
      next
    }
    state0 <- rmtif_state_at(control, midpoint)
    state1 <- rmtif_state_at(treatment, midpoint)
    if (is.na(state0) || is.na(state1)) {
      next
    }
    weight <- 1
    if (weighted) {
      argument0 <- min(midpoint, control$death_time)
      argument1 <- min(midpoint, treatment$death_time)
      survival0 <- .rmtif_km_left(km0, argument0)
      survival1 <- .rmtif_km_left(km1, argument1)
      if (survival0 <= 0 || survival1 <= 0 ||
          survival0 < minimum_survival || survival1 < minimum_survival) {
        stop("Censoring survival fell below minimum_survival; reduce tau or " ,
             "justify a different positivity threshold.", call. = FALSE)
      }
      weight <- 1 / (survival0 * survival1)
    }
    if (state1 < state0) {
      win <- win + width * weight
    } else if (state1 > state0) {
      loss <- loss + width * weight
    }
  }
  c(win = win, loss = loss, net = win - loss)
}

#' Compute pairwise restricted win, loss, and net times
#'
#' This is the exact interval-by-interval calculation for two piecewise-constant
#' paths. With `censoring = "none"`, both paths must be complete through
#' `tau`. With `censoring = "arm_km"`, the grid also contains every censoring
#' time from `reference_subjects`, and the integrand uses the terminal-aware,
#' arm-specific inverse censoring weight described in the revised manuscript.
#'
#' @param control One arm-0 `rmtif_subject`.
#' @param treatment One arm-1 `rmtif_subject`.
#' @param tau Restriction horizon.
#' @param censoring Either `"none"` or `"arm_km"`.
#' @param reference_subjects Full sample used to estimate censoring curves.
#' @param minimum_survival Positivity threshold for each censoring survival.
#' @return Named vector with `win`, `loss`, and `net` in time units.
#' @export
rmtif_pair_times <- function(control, treatment, tau,
                             censoring = c("none", "arm_km"),
                             reference_subjects = NULL,
                             minimum_survival = 0) {
  censoring <- match.arg(censoring)
  tau <- .rmtif_scalar_numeric(tau, "tau", lower = 0)
  minimum_survival <- .rmtif_minimum_survival(minimum_survival)
  if (!inherits(control, "rmtif_subject") ||
      !inherits(treatment, "rmtif_subject") ||
      tau > min(control$tau, treatment$tau)) {
    stop("Both subjects must be observed rmtif_subject paths with horizon at least tau.",
         call. = FALSE)
  }
  if (censoring == "none") {
    if (!.rmtif_available(control, tau) ||
        !.rmtif_available(treatment, tau)) {
      stop("censoring = 'none' requires both complete paths through tau.",
           call. = FALSE)
    }
    return(.rmtif_pair_integral(control, treatment, tau))
  }
  if (is.null(reference_subjects)) {
    stop("reference_subjects is required for arm_km weighting.", call. = FALSE)
  }
  split <- .rmtif_validate_subjects(reference_subjects, tau)
  km0 <- .rmtif_censor_km(split$control, tau)
  km1 <- .rmtif_censor_km(split$treatment, tau)
  support0 <- .rmtif_km_left(km0, tau)
  support1 <- .rmtif_km_left(km1, tau)
  if (support0 <= 0 || support1 <= 0 ||
      support0 < minimum_survival || support1 < minimum_survival) {
    stop("Arm-specific censoring positivity fails through tau.",
         call. = FALSE)
  }
  .rmtif_pair_integral(
    control, treatment, tau, km0, km1,
    extra_times = c(km0$times, km1$times),
    minimum_survival = minimum_survival
  )
}

.rmtif_pair_matrices <- function(subjects, tau,
                                 censoring = c("none", "arm_km"),
                                 minimum_survival = 0) {
  censoring <- match.arg(censoring)
  split <- .rmtif_validate_subjects(subjects, tau, minimum_per_arm = 2L)
  if (censoring == "none") {
    incomplete <- vapply(subjects, function(x) {
      !.rmtif_available(x, tau)
    }, logical(1))
    if (any(incomplete)) {
      stop("Complete-data estimation requires every path through tau; use " ,
           "censoring = 'arm_km' for loss to follow-up.", call. = FALSE)
    }
    km0 <- km1 <- NULL
    extra_times <- numeric()
  } else {
    km0 <- .rmtif_censor_km(split$control, tau)
    km1 <- .rmtif_censor_km(split$treatment, tau)
    support0 <- .rmtif_km_left(km0, tau)
    support1 <- .rmtif_km_left(km1, tau)
    if (support0 <= 0 || support1 <= 0 ||
        support0 < minimum_survival || support1 < minimum_survival) {
      stop("Arm-specific censoring positivity fails through tau.",
           call. = FALSE)
    }
    extra_times <- c(km0$times, km1$times)
  }
  win <- loss <- matrix(0, nrow = split$n0, ncol = split$n1)
  grid <- sort(unique(c(
    0, tau, extra_times,
    unlist(lapply(subjects, `[[`, "times"), use.names = FALSE),
    vapply(subjects, function(x) {
      if (isTRUE(x$censor_event)) x$censor_time else Inf
    }, numeric(1)),
    vapply(subjects, `[[`, numeric(1), "death_time")
  )))
  grid <- grid[is.finite(grid) & grid >= 0 & grid <= tau]
  if (censoring == "none") {
    # Evaluate each complete subject once on the full grid, rather than once
    # per interval. Interval widths, strict comparisons, and accumulation
    # order are unchanged, so this is the same exact pairwise integration.
    midpoints <- (utils::head(grid, -1L) + utils::tail(grid, -1L)) / 2
    state_matrix <- function(group) {
      value <- vapply(group, rmtif_state_at, numeric(length(midpoints)),
                      time = midpoints)
      t(matrix(value, nrow = length(midpoints), ncol = length(group)))
    }
    states0 <- state_matrix(split$control)
    states1 <- state_matrix(split$treatment)
    for (k in seq_len(length(grid) - 1L)) {
      width <- grid[k + 1L] - grid[k]
      win <- win + width * outer(states0[, k], states1[, k],
                                function(x, y) y < x)
      loss <- loss + width * outer(states0[, k], states1[, k],
                                  function(x, y) y > x)
    }
    return(list(win = win, loss = loss, net = win - loss, split = split,
                km0 = NULL, km1 = NULL, censoring = censoring))
  }
  for (k in seq_len(length(grid) - 1L)) {
    midpoint <- (grid[k] + grid[k + 1L]) / 2
    width <- grid[k + 1L] - grid[k]
    state0 <- vapply(split$control, rmtif_state_at, numeric(1),
                     time = midpoint)
    state1 <- vapply(split$treatment, rmtif_state_at, numeric(1),
                     time = midpoint)
    available0 <- is.finite(state0)
    available1 <- is.finite(state1)
    weight0 <- as.numeric(available0)
    weight1 <- as.numeric(available1)
    if (censoring == "arm_km") {
      argument0 <- vapply(split$control, function(x) {
        min(midpoint, x$death_time)
      }, numeric(1))
      argument1 <- vapply(split$treatment, function(x) {
        min(midpoint, x$death_time)
      }, numeric(1))
      survival0 <- .rmtif_km_left(km0, argument0)
      survival1 <- .rmtif_km_left(km1, argument1)
      if (any(survival0[available0] <= 0) ||
          any(survival1[available1] <= 0) ||
          any(survival0[available0] < minimum_survival) ||
          any(survival1[available1] < minimum_survival)) {
        stop("Censoring survival fell below minimum_survival; reduce tau or " ,
             "justify a different positivity threshold.", call. = FALSE)
      }
      weight0[available0] <- 1 / survival0[available0]
      weight1[available1] <- 1 / survival1[available1]
    }
    state0[!available0] <- 0
    state1[!available1] <- 0
    pair_weight <- width * outer(weight0, weight1)
    win <- win + pair_weight * outer(state0, state1, function(x, y) y < x)
    loss <- loss + pair_weight * outer(state0, state1, function(x, y) y > x)
  }
  list(win = win, loss = loss, net = win - loss, split = split,
       km0 = km0, km1 = km1, censoring = censoring)
}

# Exact state-mass aggregation for the IPCW mean --------------------------------
#
# The marginal IPCW estimator is the average of all pair contributions.  On a
# fixed open integration interval, that average can be evaluated from the
# arm-specific weighted state masses rather than by materializing every
# n0-by-n1 pair.  This is algebraically identical to the strict pairwise
# comparison used in .rmtif_pair_matrices(): for ordered states r and q,
# treatment wins when q < r and loses when q > r.  It is used only where a
# scalar mean is needed (not for the fixed-weight projection diagnostic).

.rmtif_weighted_state_pair_totals <- function(state0, weight0, state1, weight1) {
  keep0 <- is.finite(state0) & is.finite(weight0) & weight0 > 0
  keep1 <- is.finite(state1) & is.finite(weight1) & weight1 > 0
  if (!any(keep0) || !any(keep1)) {
    return(c(win = 0, loss = 0, net = 0))
  }
  ordered_states <- sort(unique(c(state0[keep0], state1[keep1])))
  mass0 <- vapply(ordered_states, function(value) {
    sum(weight0[keep0 & state0 == value])
  }, numeric(1))
  mass1 <- vapply(ordered_states, function(value) {
    sum(weight1[keep1 & state1 == value])
  }, numeric(1))
  below_treatment <- c(0, utils::head(cumsum(mass1), -1L))
  above_treatment <- sum(mass1) - cumsum(mass1)
  win <- sum(mass0 * below_treatment)
  loss <- sum(mass0 * above_treatment)
  c(win = win, loss = loss, net = win - loss)
}

.rmtif_ipcw_mean_fast <- function(subjects, tau, minimum_survival = 0) {
  tau <- .rmtif_scalar_numeric(tau, "tau", lower = 0)
  minimum_survival <- .rmtif_minimum_survival(minimum_survival)
  split <- .rmtif_validate_subjects(subjects, tau)
  km0 <- .rmtif_censor_km(split$control, tau)
  km1 <- .rmtif_censor_km(split$treatment, tau)
  support0 <- .rmtif_km_left(km0, tau)
  support1 <- .rmtif_km_left(km1, tau)
  if (support0 <= 0 || support1 <= 0 ||
      support0 < minimum_survival || support1 < minimum_survival) {
    stop("Arm-specific censoring positivity fails through tau.",
         call. = FALSE)
  }

  grid <- sort(unique(c(
    0, tau, km0$times, km1$times,
    unlist(lapply(subjects, `[[`, "times"), use.names = FALSE),
    vapply(subjects, function(x) {
      if (isTRUE(x$censor_event)) x$censor_time else Inf
    }, numeric(1)),
    vapply(subjects, `[[`, numeric(1), "death_time")
  )))
  grid <- grid[is.finite(grid) & grid >= 0 & grid <= tau]
  total <- c(win = 0, loss = 0, net = 0)
  if (length(grid) < 2L) {
    return(total)
  }

  for (interval in seq_len(length(grid) - 1L)) {
    midpoint <- (grid[interval] + grid[interval + 1L]) / 2
    width <- grid[interval + 1L] - grid[interval]
    state0 <- vapply(split$control, rmtif_state_at, numeric(1),
                     time = midpoint)
    state1 <- vapply(split$treatment, rmtif_state_at, numeric(1),
                     time = midpoint)
    available0 <- is.finite(state0)
    available1 <- is.finite(state1)
    if (!any(available0) || !any(available1)) {
      next
    }
    argument0 <- vapply(split$control, function(x) {
      min(midpoint, x$death_time)
    }, numeric(1))
    argument1 <- vapply(split$treatment, function(x) {
      min(midpoint, x$death_time)
    }, numeric(1))
    survival0 <- .rmtif_km_left(km0, argument0)
    survival1 <- .rmtif_km_left(km1, argument1)
    if (any(survival0[available0] <= 0) ||
        any(survival1[available1] <= 0) ||
        any(survival0[available0] < minimum_survival) ||
        any(survival1[available1] < minimum_survival)) {
      stop("Censoring survival fell below minimum_survival; reduce tau or ",
           "justify a different positivity threshold.", call. = FALSE)
    }
    weight0 <- numeric(length(state0))
    weight1 <- numeric(length(state1))
    weight0[available0] <- 1 / survival0[available0]
    weight1[available1] <- 1 / survival1[available1]
    total <- total + width * .rmtif_weighted_state_pair_totals(
      state0, weight0, state1, weight1
    )
  }
  total / (split$n0 * split$n1)
}

# Bootstrap cache for the scalar IPCW mean -------------------------------------
#
# A bootstrap sample only changes subject multiplicities and the two fitted KM
# curves.  Its state, death, and censoring times are a subset of the original
# union grid.  Consequently, a grid based on the original records is an exact
# refinement of every resample's grid.  Caching the states on that grid and
# applying multinomial resampling counts avoids rebuilding trajectories or
# re-evaluating paths for every bootstrap replicate.

.rmtif_ipcw_bootstrap_cache <- function(subjects, tau) {
  tau <- .rmtif_scalar_numeric(tau, "tau", lower = 0)
  split <- .rmtif_validate_subjects(subjects, tau, minimum_per_arm = 2L)
  grid <- sort(unique(c(
    0, tau,
    unlist(lapply(subjects, `[[`, "times"), use.names = FALSE),
    vapply(subjects, function(x) {
      if (isTRUE(x$censor_event)) x$censor_time else Inf
    }, numeric(1)),
    vapply(subjects, `[[`, numeric(1), "death_time")
  )))
  grid <- grid[is.finite(grid) & grid >= 0 & grid <= tau]
  midpoint <- (utils::head(grid, -1L) + utils::tail(grid, -1L)) / 2
  state_matrix <- function(group) {
    if (!length(midpoint)) {
      return(matrix(numeric(), nrow = length(group), ncol = 0L))
    }
    values <- vapply(group, rmtif_state_at, numeric(length(midpoint)),
                     time = midpoint)
    matrix(values, nrow = length(group), ncol = length(midpoint),
           byrow = TRUE)
  }
  death0 <- vapply(split$control, `[[`, numeric(1), "death_time")
  death1 <- vapply(split$treatment, `[[`, numeric(1), "death_time")
  argument_matrix <- function(death) {
    matrix(pmin(rep(midpoint, each = length(death)),
                rep(death, length(midpoint))), nrow = length(death))
  }
  list(
    split = split,
    tau = tau,
    midpoint = midpoint,
    width = diff(grid),
    state0 = state_matrix(split$control),
    state1 = state_matrix(split$treatment),
    death0 = death0,
    death1 = death1,
    argument0 = argument_matrix(death0),
    argument1 = argument_matrix(death1),
    ordered_states = sort(unique(c(
      unlist(lapply(split$control, `[[`, "states"), use.names = FALSE),
      unlist(lapply(split$treatment, `[[`, "states"), use.names = FALSE)
    )))
  )
}

.rmtif_validate_bootstrap_indices <- function(index, n, arm_label) {
  if (!is.numeric(index) || !length(index) || anyNA(index) ||
      any(!is.finite(index)) || any(index != floor(index)) ||
      any(index < 1L) || any(index > n)) {
    stop(
      arm_label, " bootstrap indices must be finite whole numbers from 1 through ",
      n, ".", call. = FALSE
    )
  }
  as.integer(index)
}

.rmtif_ipcw_bootstrap_mean_from_cache <- function(cache, control_index,
                                                   treatment_index,
                                                   minimum_survival = 0) {
  minimum_survival <- .rmtif_minimum_survival(minimum_survival)
  control_index <- .rmtif_validate_bootstrap_indices(
    control_index, cache$split$n0, "Control"
  )
  treatment_index <- .rmtif_validate_bootstrap_indices(
    treatment_index, cache$split$n1, "Treatment"
  )
  n0 <- length(control_index)
  n1 <- length(treatment_index)
  if (!n0 || !n1) {
    stop("Each bootstrap arm must contain at least one sampled subject.",
         call. = FALSE)
  }
  counts0 <- tabulate(control_index, nbins = cache$split$n0)
  counts1 <- tabulate(treatment_index, nbins = cache$split$n1)
  km0 <- .rmtif_censor_km(cache$split$control[control_index], cache$tau)
  km1 <- .rmtif_censor_km(cache$split$treatment[treatment_index], cache$tau)
  support0 <- .rmtif_km_left(km0, cache$tau)
  support1 <- .rmtif_km_left(km1, cache$tau)
  if (support0 <= 0 || support1 <= 0 ||
      support0 < minimum_survival || support1 < minimum_survival) {
    stop("Arm-specific censoring positivity fails through tau.",
         call. = FALSE)
  }
  m <- length(cache$midpoint)
  number_of_states <- length(cache$ordered_states)
  if (!m || !number_of_states) {
    return(c(win = 0, loss = 0, net = 0))
  }
  # Each column is one exact union-grid interval. All KM arguments are
  # predictable left limits at min(midpoint, observed death), as in the
  # pairwise formula. The KM is refitted on the resampled subjects above.
  survival0 <- matrix(.rmtif_km_left(km0, cache$argument0),
                      nrow = cache$split$n0)
  survival1 <- matrix(.rmtif_km_left(km1, cache$argument1),
                      nrow = cache$split$n1)
  weight0 <- counts0 / survival0
  weight1 <- counts1 / survival1
  mass_by_state <- function(state, weight) {
    t(matrix(vapply(cache$ordered_states, function(level) {
      colSums(weight * (state == level), na.rm = TRUE)
    }, numeric(m)), nrow = m, ncol = number_of_states))
  }
  mass0 <- mass_by_state(cache$state0, weight0)
  mass1 <- mass_by_state(cache$state1, weight1)
  cumulative1 <- matrix(apply(mass1, 2L, cumsum),
                        nrow = number_of_states, ncol = m)
  below1 <- cumulative1 - mass1
  above1 <- matrix(rep(colSums(mass1), each = number_of_states),
                   nrow = number_of_states) - cumulative1
  win <- sum(cache$width * colSums(mass0 * below1)) / (n0 * n1)
  loss <- sum(cache$width * colSums(mass0 * above1)) / (n0 * n1)
  c(win = win, loss = loss, net = win - loss)
}

.rmtif_projection_from_components <- function(estimate, row_projection,
                                               column_projection,
                                               conf_level) {
  conf_level <- .rmtif_conf_level(conf_level)
  n0 <- length(row_projection)
  n1 <- length(column_projection)
  if (n0 < 2L || n1 < 2L) {
    stop("Projection-variance inference requires at least two subjects per arm.",
         call. = FALSE)
  }
  variance <- stats::var(row_projection) / n0 +
    stats::var(column_projection) / n1
  se <- sqrt(max(0, variance))
  critical <- stats::qnorm(1 - (1 - conf_level) / 2)
  statistic <- if (se > 0) estimate / se else NA_real_
  p_value <- if (is.finite(statistic)) {
    2 * stats::pnorm(-abs(statistic))
  } else {
    NA_real_
  }
  list(
    estimate = estimate,
    se = se,
    conf_int = c(lower = estimate - critical * se,
                 upper = estimate + critical * se),
    statistic = statistic,
    p_value = p_value,
    row_projection = row_projection,
    column_projection = column_projection
  )
}

.rmtif_projection_summary <- function(kernel, conf_level) {
  .rmtif_projection_from_components(
    estimate = mean(kernel),
    row_projection = rowMeans(kernel),
    column_projection = colMeans(kernel),
    conf_level = conf_level
  )
}

#' Compute two-sample projection inference from pairwise net times
#'
#' The estimate is the mean of a control-by-treatment net-time matrix. Its
#' standard error combines the empirical variance of row means divided by the
#' control sample size with the variance of column means divided by the
#' treatment sample size. A normal confidence interval and two-sided test of
#' zero contrast use the supplied confidence level.
#'
#' @param pairwise_net_times Finite numeric matrix with control subjects in
#'   rows and treatment subjects in columns.
#' @param conf_level Confidence level for the normal interval.
#' @return List with `estimate`, `se`, `conf_int`, `statistic`, `p_value`,
#'   `row_projection`, and `column_projection`.
#' @export
rmtif_projection_inference <- function(pairwise_net_times, conf_level) {
  pairwise_net_times <- as.matrix(pairwise_net_times)
  if (!is.numeric(pairwise_net_times) || anyNA(pairwise_net_times) ||
      any(!is.finite(pairwise_net_times)) ||
      nrow(pairwise_net_times) < 2L || ncol(pairwise_net_times) < 2L) {
    stop("pairwise_net_times must be a finite numeric matrix with at least ",
         "two subjects in each arm.", call. = FALSE)
  }
  .rmtif_projection_summary(pairwise_net_times, conf_level)
}

#' Estimate complete-data rmt-IF with projection inference
#'
#' Average pairwise treatment-minus-control net times over every control and
#' treatment subject pair. Subject-level row and column projections yield the
#' two-sample U-statistic standard error. The function reports a normal
#' confidence interval and the z statistic and two-sided p-value for a zero
#' contrast; no independent-increment assumption across time is used.
#'
#' @param subjects List of `rmtif_subject` objects with complete paths to `tau`.
#' @param tau Restriction horizon.
#' @param conf_level Confidence level.
#' @return An `rmtif_fit` list containing `estimate`, projection `se`,
#'   `conf_int`, zero-contrast `statistic` and `p_value`, mean pairwise `win`
#'   and `loss` times, and row and column projections.
#' @export
rmtif_complete <- function(subjects, tau, conf_level) {
  conf_level <- .rmtif_conf_level(conf_level)
  pair <- .rmtif_pair_matrices(subjects, tau, censoring = "none")
  result <- rmtif_projection_inference(pair$net, conf_level)
  result$method <- "Complete-data two-sample projection"
  result$n0 <- pair$split$n0
  result$n1 <- pair$split$n1
  result$win <- mean(pair$win)
  result$loss <- mean(pair$loss)
  class(result) <- "rmtif_fit"
  result
}

.rmtif_ipcw_point <- function(subjects, tau, minimum_survival, conf_level) {
  pair <- .rmtif_pair_matrices(subjects, tau, censoring = "arm_km",
                               minimum_survival = minimum_survival)
  projection <- .rmtif_projection_summary(pair$net, conf_level)
  list(
    estimate = mean(pair$net),
    win = mean(pair$win),
    loss = mean(pair$loss),
    projection_diagnostic = projection,
    pair = pair
  )
}

#' Check the exact IPCW state-mass implementation against dense pair averaging
#'
#' This audit helper evaluates one observed data set by both algebraically
#' identical routes.  It is useful for regression testing the faster bootstrap
#' calculation; it does not fit a bootstrap or change the estimator.
#'
#' @param subjects List of observed `rmtif_subject` objects.
#' @param tau Restriction horizon.
#' @param minimum_survival Positivity threshold for each arm-specific curve.
#' @param tolerance Absolute numerical comparison tolerance.
#' @return One-row data frame with dense, direct state-mass, and cached
#'   state-mass IPCW means, including a deterministic duplicate-subject
#'   bootstrap resample.
#' @export
rmtif_ipcw_equivalence_check <- function(subjects, tau,
                                         minimum_survival = 0,
                                         tolerance = 1e-10) {
  tolerance <- .rmtif_scalar_numeric(tolerance, "tolerance", lower = 0)
  dense <- mean(.rmtif_pair_matrices(
    subjects, tau, censoring = "arm_km",
    minimum_survival = minimum_survival
  )$net)
  state_mass <- .rmtif_ipcw_mean_fast(
    subjects, tau, minimum_survival
  )["net"]
  cache <- .rmtif_ipcw_bootstrap_cache(subjects, tau)
  cached_state_mass <- .rmtif_ipcw_bootstrap_mean_from_cache(
    cache,
    control_index = seq_len(cache$split$n0),
    treatment_index = seq_len(cache$split$n1),
    minimum_survival = minimum_survival
  )["net"]
  duplicate_index <- function(n) {
    c(1L, 1L, if (n > 2L) seq.int(3L, n) else integer())
  }
  bootstrap_control_index <- duplicate_index(cache$split$n0)
  bootstrap_treatment_index <- duplicate_index(cache$split$n1)
  dense_bootstrap <- mean(.rmtif_pair_matrices(
    c(
      cache$split$control[bootstrap_control_index],
      cache$split$treatment[bootstrap_treatment_index]
    ),
    tau, censoring = "arm_km", minimum_survival = minimum_survival
  )$net)
  cached_bootstrap <- .rmtif_ipcw_bootstrap_mean_from_cache(
    cache, bootstrap_control_index, bootstrap_treatment_index,
    minimum_survival
  )["net"]
  difference <- as.numeric(dense - state_mass)
  cached_difference <- as.numeric(dense - cached_state_mass)
  bootstrap_difference <- as.numeric(dense_bootstrap - cached_bootstrap)
  if (!is.finite(difference) || !is.finite(cached_difference) ||
      !is.finite(bootstrap_difference) || abs(difference) > tolerance ||
      abs(cached_difference) > tolerance || abs(bootstrap_difference) > tolerance) {
    stop(
      "Dense-pair and state-mass IPCW means disagree (direct difference = ",
      format(difference, digits = 16), "; cached difference = ",
      format(cached_difference, digits = 16), "; bootstrap difference = ",
      format(bootstrap_difference, digits = 16), ").",
      call. = FALSE
    )
  }
  data.frame(
    dense_pair_mean = dense,
    state_mass_mean = as.numeric(state_mass),
    cached_state_mass_mean = as.numeric(cached_state_mass),
    bootstrap_dense_pair_mean = dense_bootstrap,
    bootstrap_cached_state_mass_mean = as.numeric(cached_bootstrap),
    absolute_difference = abs(difference),
    cached_absolute_difference = abs(cached_difference),
    bootstrap_cached_absolute_difference = abs(bootstrap_difference),
    row.names = NULL
  )
}

#' Audit arm-specific censoring support for IPCW
#'
#' The aim is to check the empirical positivity condition before estimation.
#' Survival is evaluated as the predictable left limit at `tau`; therefore a
#' Kaplan--Meier curve that first reaches zero exactly at `tau` still supports
#' integration through that endpoint, whereas a zero before `tau` does not.
#'
#' @param subjects List of `rmtif_subject` objects.
#' @param tau Candidate restriction horizon.
#' @return One row per arm with the left-limit survival at `tau`, the first
#'   zero observed by `tau`, maximum empirically verified horizon, and a
#'   support flag. Absence of a zero by `tau` is not extrapolated beyond `tau`.
#' @export
rmtif_censoring_support <- function(subjects, tau) {
  split <- .rmtif_validate_subjects(subjects, tau)
  summaries <- lapply(list(split$control, split$treatment), function(group) {
    km <- .rmtif_censor_km(group, tau)
    zero_index <- which(km$survival <= 0)
    first_zero <- if (length(zero_index)) km$times[zero_index[1L]] else Inf
    data.frame(
      survival_left_at_tau = .rmtif_km_left(km, tau),
      first_zero_time = first_zero,
      maximum_verified_horizon = if (is.finite(first_zero)) first_zero else tau,
      supported = tau <= first_zero,
      row.names = NULL
    )
  })
  result <- do.call(rbind, summaries)
  result$arm <- 0:1
  result <- result[, c("arm", setdiff(names(result), "arm"))]
  rownames(result) <- NULL
  result
}

.rmtif_bootstrap_samples <- function(split, bootstrap, statistic) {
  values <- rep(NA_real_, bootstrap)
  for (b in seq_len(bootstrap)) {
    sample0 <- split$control[sample.int(split$n0, split$n0, replace = TRUE)]
    sample1 <- split$treatment[sample.int(split$n1, split$n1, replace = TRUE)]
    values[b] <- tryCatch(statistic(c(sample0, sample1)),
                          error = function(e) NA_real_)
  }
  values
}

#' Estimate terminal-aware arm-specific IPCW rmt-IF
#'
#' The aim is unbiased marginal rmt-IF estimation under unconditional
#' independent loss to follow-up within treatment arm.  The estimator uses two
#' arm-specific censoring Kaplan-Meier curves, carries observed death forward,
#' and integrates over the union of state, death, and censoring times.  Primary
#' inference is a treatment-stratified bootstrap that refits both censoring
#' curves in every replicate.  The fixed-weight projection is returned only as
#' a diagnostic because it omits censoring-curve estimation uncertainty.
#'
#' @param subjects List of observed `rmtif_subject` objects.
#' @param tau Restriction horizon.
#' @param bootstrap Number of stratified bootstrap replicates.
#' @param conf_level Confidence level.
#' @param bootstrap_success_fraction Minimum fraction of finite bootstrap
#'   estimates needed to report bootstrap inference.
#' @param seed Optional reproducibility seed.
#' @param minimum_survival Positivity threshold for each arm-specific curve.
#' @param projection_diagnostic Whether to form the dense original-sample pair
#'   matrix for the fixed-weight projection diagnostic.  Set to `FALSE` in a
#'   Monte Carlo study when only the IPCW point estimate and refitted-bootstrap
#'   inference are required.
#' @param allow_insufficient_bootstrap If `TRUE`, return a diagnostic object
#'   when too few bootstrap resamples retain support.  Its bootstrap
#'   standard error and intervals are left undefined; this option is intended
#'   only to audit a boundary horizon, not to report primary inference.
#' @return An `rmtif_fit` list with bootstrap inference when requested.
#' @export
rmtif_ipcw <- function(subjects, tau, bootstrap, conf_level,
                       bootstrap_success_fraction, seed = NULL,
                       minimum_survival = 0,
                       projection_diagnostic = TRUE,
                       allow_insufficient_bootstrap = FALSE) {
  tau <- .rmtif_scalar_numeric(tau, "tau", lower = 0)
  conf_level <- .rmtif_conf_level(conf_level)
  bootstrap_success_fraction <- .rmtif_bootstrap_success_fraction(
    bootstrap_success_fraction
  )
  minimum_survival <- .rmtif_minimum_survival(minimum_survival)
  bootstrap <- .rmtif_scalar_count(bootstrap, "bootstrap", lower = 0)
  if (bootstrap == 1L) {
    stop("bootstrap must be 0 or at least 2.", call. = FALSE)
  }
  projection_diagnostic <- .rmtif_scalar_flag(
    projection_diagnostic, "projection_diagnostic"
  )
  allow_insufficient_bootstrap <- .rmtif_scalar_flag(
    allow_insufficient_bootstrap, "allow_insufficient_bootstrap"
  )
  split <- .rmtif_validate_subjects(subjects, tau, minimum_per_arm = 2L)
  cache <- if (bootstrap >= 2L || !projection_diagnostic) {
    .rmtif_ipcw_bootstrap_cache(subjects, tau)
  } else NULL
  point <- if (projection_diagnostic) {
    .rmtif_ipcw_point(subjects, tau, minimum_survival, conf_level)
  } else {
    state_mass <- .rmtif_ipcw_bootstrap_mean_from_cache(
      cache, seq_len(split$n0), seq_len(split$n1), minimum_survival
    )
    list(
      estimate = unname(state_mass["net"]),
      win = unname(state_mass["win"]),
      loss = unname(state_mass["loss"]),
      projection_diagnostic = NULL
    )
  }
  result <- list(
    estimate = point$estimate,
    win = point$win,
    loss = point$loss,
    n0 = split$n0,
    n1 = split$n1,
    projection_diagnostic = point$projection_diagnostic,
    bootstrap_estimates = numeric(),
    bootstrap_attempted = bootstrap,
    bootstrap_success = 0L,
    bootstrap_supported = NA,
    bootstrap_message = "Bootstrap not requested.",
    bootstrap_failures = character(),
    se = NA_real_,
    conf_int = c(lower = NA_real_, upper = NA_real_),
    percentile_conf_int = c(lower = NA_real_, upper = NA_real_),
    statistic = NA_real_,
    p_value = NA_real_,
    method = "Arm-specific terminal-aware IPCW"
  )
  if (bootstrap >= 2L) {
    bootstrap_result <- .rmtif_with_seed(seed, {
      result_values <- rep(NA_real_, bootstrap)
      result_errors <- rep("", bootstrap)
      for (replicate in seq_len(bootstrap)) {
        control_index <- sample.int(split$n0, split$n0, replace = TRUE)
        treatment_index <- sample.int(split$n1, split$n1, replace = TRUE)
        outcome <- tryCatch(
          list(value = as.numeric(.rmtif_ipcw_bootstrap_mean_from_cache(
            cache, control_index, treatment_index, minimum_survival
          )["net"]), error = ""),
          error = function(e) list(value = NA_real_, error = conditionMessage(e))
        )
        result_values[replicate] <- outcome$value
        result_errors[replicate] <- outcome$error
      }
      list(values = result_values, errors = result_errors)
    })
    values <- bootstrap_result$values
    failure_messages <- unique(bootstrap_result$errors[
      nzchar(bootstrap_result$errors)
    ])
    good <- values[is.finite(values)]
    result$bootstrap_estimates <- values
    result$bootstrap_success <- length(good)
    result$bootstrap_failures <- failure_messages
    if (length(good) < max(2L, ceiling(
      bootstrap_success_fraction * bootstrap
    ))) {
      result$bootstrap_supported <- FALSE
      support_failures <- length(failure_messages) > 0L && all(grepl(
        "positivity|minimum_survival", failure_messages, ignore.case = TRUE
      ))
      result$bootstrap_message <- if (support_failures) {
        paste0(
          "Only ", length(good), "/", bootstrap,
          " IPCW bootstrap resamples retained censoring support; primary ",
          "bootstrap inference is not reported at this horizon."
        )
      } else {
        paste0(
          "Only ", length(good), "/", bootstrap,
          " IPCW bootstrap resamples produced finite estimates; primary ",
          "bootstrap inference is not reported. Inspect bootstrap_failures."
        )
      }
      if (!isTRUE(allow_insufficient_bootstrap)) {
        stop(result$bootstrap_message, call. = FALSE)
      }
      class(result) <- "rmtif_fit"
      return(result)
    }
    result$bootstrap_supported <- TRUE
    result$bootstrap_message <- paste0(
      length(good), "/", bootstrap,
      " IPCW bootstrap resamples produced finite estimates with refitted ",
      "censoring curves."
    )
    result$se <- stats::sd(good)
    critical <- stats::qnorm(1 - (1 - conf_level) / 2)
    result$conf_int <- c(lower = result$estimate - critical * result$se,
                         upper = result$estimate + critical * result$se)
    result$percentile_conf_int <- stats::quantile(
      good, probs = c((1 - conf_level) / 2, 1 - (1 - conf_level) / 2),
      names = FALSE, type = 6
    )
    names(result$percentile_conf_int) <- c("lower", "upper")
    if (is.finite(result$se) && result$se > 0) {
      result$statistic <- result$estimate / result$se
      result$p_value <- 2 * stats::pnorm(-abs(result$statistic))
    }
  }
  class(result) <- "rmtif_fit"
  result
}

#' Estimate a complete-data rmt-IF curve
#'
#' The function generates pointwise estimates and projection intervals for
#' complete ordered-state paths.
#'
#' @param subjects Complete observed paths.
#' @param times Increasing restriction horizons.
#' @param conf_level Pointwise confidence level.
#' @return Data frame with estimate, standard error, interval, and p-value.
#' @export
rmtif_complete_curve <- function(subjects, times, conf_level) {
  conf_level <- .rmtif_conf_level(conf_level)
  times <- sort(unique(as.numeric(times)))
  if (!length(times) || anyNA(times) || any(!is.finite(times)) ||
      any(times <= 0)) {
    stop("times must contain positive finite restriction horizons.", call. = FALSE)
  }
  tau <- max(times)
  split <- .rmtif_validate_subjects(subjects, tau, minimum_per_arm = 2L)
  incomplete <- vapply(subjects, function(x) !.rmtif_available(x, tau),
                       logical(1))
  if (any(incomplete)) {
    stop("Complete-data curves require every path through max(times).",
         call. = FALSE)
  }
  grid <- sort(unique(c(
    0, times,
    unlist(lapply(subjects, `[[`, "times"), use.names = FALSE),
    vapply(subjects, `[[`, numeric(1), "death_time")
  )))
  grid <- grid[is.finite(grid) & grid >= 0 & grid <= tau]
  # The dense pair matrix is unnecessary for a curve.  On each interval the
  # exact row and column kernel sums are obtained from strict empirical-CDF
  # counts in the opposite arm.  Accumulating those subject projections is
  # algebraically identical to updating every cell of the n0-by-n1 matrix, but
  # reduces each interval from O(n0*n1) work to sorting plus linear lookups.
  row_kernel_sum <- numeric(split$n0)
  column_kernel_sum <- numeric(split$n1)
  result <- vector("list", length(times))
  midpoints <- (utils::head(grid, -1L) + utils::tail(grid, -1L)) / 2
  state_matrix <- function(group) {
    values <- vapply(group, rmtif_state_at, numeric(length(midpoints)),
                     time = midpoints)
    t(matrix(values, nrow = length(midpoints), ncol = length(group)))
  }
  states0 <- state_matrix(split$control)
  states1 <- state_matrix(split$treatment)
  targets <- match(utils::tail(grid, -1L), times)
  for (k in seq_len(length(grid) - 1L)) {
    width <- grid[k + 1L] - grid[k]
    state0 <- states0[, k]
    state1 <- states1[, k]
    sorted0 <- sort(state0)
    sorted1 <- sort(state1)
    treatment_below_control <- findInterval(
      state0, sorted1, left.open = TRUE
    )
    treatment_above_control <- split$n1 - findInterval(state0, sorted1)
    control_above_treatment <- split$n0 - findInterval(state1, sorted0)
    control_below_treatment <- findInterval(
      state1, sorted0, left.open = TRUE
    )
    row_kernel_sum <- row_kernel_sum + width * (
      treatment_below_control - treatment_above_control
    )
    column_kernel_sum <- column_kernel_sum + width * (
      control_above_treatment - control_below_treatment
    )
    target <- targets[k]
    if (!is.na(target)) {
      row_projection <- row_kernel_sum / split$n1
      column_projection <- column_kernel_sum / split$n0
      result[[target]] <- .rmtif_projection_from_components(
        estimate = mean(row_projection),
        row_projection = row_projection,
        column_projection = column_projection,
        conf_level = conf_level
      )
    }
  }
  data.frame(
    time = times,
    estimate = vapply(result, `[[`, numeric(1), "estimate"),
    se = vapply(result, `[[`, numeric(1), "se"),
    lower = vapply(result, function(x) x$conf_int["lower"], numeric(1)),
    upper = vapply(result, function(x) x$conf_int["upper"], numeric(1)),
    p_value = vapply(result, `[[`, numeric(1), "p_value"),
    row.names = NULL
  )
}

#' Evaluate the exact finite-state projection curve by state counts
#'
#' This is the fast implementation of the complete-path pairwise estimator.
#' For each interval it obtains row and column projections from strict
#' empirical state counts in the opposite arm, algebraically matching the
#' dense pair matrix without assuming independent increments across time.
#'
#' @param subjects Complete observed ordered-state paths.
#' @param times Increasing positive restriction horizons.
#' @param conf_level Pointwise normal confidence level.
#' @return Data frame with horizon, rmt-IF estimate, projection SE, interval
#'   endpoints, and two-sided zero-contrast p-value.
#' @export
rmtif_complete_state_count_curve <- function(subjects, times, conf_level) {
  rmtif_complete_curve(subjects, times, conf_level)
}

#' Summarize Monte Carlo estimator and interval performance
#'
#' For each requested design cell, calculate the mean estimate, bias,
#' empirical standard deviation, mean reported standard error, interval
#' coverage of the supplied truth, two-sided zero-test rejection fraction,
#' and successful-fit count. Additional trial-level columns can be averaged
#' over all generated trials, as with observed loss-to-follow-up fractions.
#'
#' @param replicates Data frame with one row per simulated trial.
#' @param group_columns Character names of the columns defining design cells.
#' @param estimate_column Name of the numeric estimate column.
#' @param se_column Name of the numeric estimated-standard-error column.
#' @param truth_column Name of the numeric population-truth column.
#' @param lower_column Name of the lower confidence-limit column.
#' @param upper_column Name of the upper confidence-limit column.
#' @param p_value_column Name of the two-sided zero-test p-value column.
#' @param test_level Rejection threshold for the p-value.
#' @param extra_mean_columns Optional trial-level numeric columns to average.
#' @return One row per design cell with truth, mean estimate, bias, empirical
#'   SD, mean estimated SE, coverage, rejection, successful-fit count, total
#'   trial count, and means of requested extra columns.
#' @export
rmtif_monte_carlo_summary <- function(
    replicates, group_columns, estimate_column, se_column, truth_column,
    lower_column, upper_column, p_value_column, test_level,
    extra_mean_columns = character()) {
  if (!is.data.frame(replicates) || !nrow(replicates)) {
    stop("replicates must be a nonempty data frame.", call. = FALSE)
  }
  column_names <- c(group_columns, estimate_column, se_column, truth_column,
                    lower_column, upper_column, p_value_column,
                    extra_mean_columns)
  if (!all(column_names %in% names(replicates))) {
    stop("Every requested column must be present in replicates.",
         call. = FALSE)
  }
  test_level <- .rmtif_conf_level(test_level)
  key <- if (length(group_columns)) {
    do.call(paste, c(replicates[group_columns], sep = "\r"))
  } else {
    rep("all", nrow(replicates))
  }
  groups <- split(seq_len(nrow(replicates)), factor(key, levels = unique(key)))
  out <- lapply(groups, function(index) {
    x <- replicates[index, , drop = FALSE]
    truth <- x[[truth_column]][1L]
    if (!is.numeric(x[[truth_column]]) || !is.finite(truth) ||
        any(!is.finite(x[[truth_column]])) ||
        any(x[[truth_column]] != truth)) {
      stop("Each design cell must have one finite numeric truth.",
           call. = FALSE)
    }
    estimate <- x[[estimate_column]]
    se <- x[[se_column]]
    lower <- x[[lower_column]]
    upper <- x[[upper_column]]
    p_value <- x[[p_value_column]]
    valid_estimate <- is.finite(estimate)
    valid_inference <- valid_estimate & is.finite(se) &
      is.finite(lower) & is.finite(upper) & is.finite(p_value)
    estimated <- estimate[valid_estimate]
    mean_estimate <- if (length(estimated)) mean(estimated) else NA_real_
    result <- data.frame(
      true = truth,
      mean_estimate = mean_estimate,
      bias = mean_estimate - truth,
      empirical_sd = if (length(estimated) > 1L) stats::sd(estimated)
        else NA_real_,
      mean_estimated_se = if (any(valid_inference)) mean(se[valid_inference])
        else NA_real_,
      coverage = if (any(valid_inference)) {
        mean(lower[valid_inference] <= truth &
               upper[valid_inference] >= truth)
      } else NA_real_,
      rejection = if (any(valid_inference)) {
        mean(p_value[valid_inference] < test_level)
      } else NA_real_,
      successful = sum(valid_inference),
      replications = nrow(x)
    )
    if (length(group_columns)) {
      result <- cbind(x[1L, group_columns, drop = FALSE], result)
    }
    for (column in extra_mean_columns) {
      if (!is.numeric(x[[column]])) {
        stop("extra_mean_columns must be numeric.", call. = FALSE)
      }
      result[[paste0("mean_", column)]] <- mean(x[[column]])
    }
    result
  })
  result <- do.call(rbind, out)
  rownames(result) <- NULL
  result
}

# -----------------------------------------------------------------------------
# Pairwise regression with a treatment intercept and covariate differences
# -----------------------------------------------------------------------------

.rmtif_regression_design <- function(pair) {
  x0 <- lapply(pair$split$control, `[[`, "covariates")
  x1 <- lapply(pair$split$treatment, `[[`, "covariates")
  dimensions <- unique(c(vapply(x0, length, integer(1)),
                         vapply(x1, length, integer(1))))
  if (length(dimensions) != 1L) {
    stop("All subjects must have the same covariate dimension.", call. = FALSE)
  }
  p <- dimensions[1L]
  matrix0 <- if (p) matrix(unlist(x0), nrow = pair$split$n0, byrow = TRUE) else
    matrix(numeric(), nrow = pair$split$n0, ncol = 0L)
  matrix1 <- if (p) matrix(unlist(x1), nrow = pair$split$n1, byrow = TRUE) else
    matrix(numeric(), nrow = pair$split$n1, ncol = 0L)
  list(
    z0 = cbind(treatment = 0, matrix0),
    z1 = cbind(treatment = 1, matrix1),
    names = c("treatment", if (p) paste0("covariate", seq_len(p)) else character())
  )
}

.rmtif_regression_components <- function(beta, win, loss, z0, z1) {
  n0 <- nrow(win)
  n1 <- ncol(win)
  eta0 <- as.numeric(z0 %*% beta)
  eta1 <- as.numeric(z1 %*% beta)
  raw_exponent <- outer(-eta0, eta1, "+")
  exponent_clipped <- any(!is.finite(raw_exponent) |
                           raw_exponent < -700 | raw_exponent > 700)
  exponent <- pmax(-700, pmin(700, raw_exponent))
  scaled_loss <- loss * exp(exponent)
  residual <- win - scaled_loss
  divisor <- n0 * n1
  row_residual <- rowSums(residual)
  column_residual <- colSums(residual)
  score <- (crossprod(z1, column_residual) -
              crossprod(z0, row_residual)) / divisor

  row_loss <- rowSums(scaled_loss)
  column_loss <- colSums(scaled_loss)
  term11 <- crossprod(z1, z1 * column_loss)
  term10 <- crossprod(z1, t(scaled_loss) %*% z0)
  term01 <- crossprod(z0, scaled_loss %*% z1)
  term00 <- crossprod(z0, z0 * row_loss)
  jacobian <- -(term11 - term10 - term01 + term00) / divisor
  list(score = as.numeric(score), jacobian = as.matrix(jacobian),
       residual = residual,
       exponent_clipped = exponent_clipped,
       maximum_abs_exponent = max(abs(raw_exponent)))
}

.rmtif_regression_from_pair <- function(pair, conf_level,
                                        beta_initial = NULL,
                                        maximum_iterations = 300L,
                                        score_tolerance = 1e-7) {
  conf_level <- .rmtif_conf_level(conf_level)
  score_tolerance <- .rmtif_scalar_numeric(
    score_tolerance, "score_tolerance", lower = 0
  )
  design <- .rmtif_regression_design(pair)
  q <- ncol(design$z0)
  if (is.null(beta_initial)) {
    beta_initial <- numeric(q)
  }
  if (length(beta_initial) != q || any(!is.finite(beta_initial))) {
    stop("beta_initial has the wrong length or non-finite values.", call. = FALSE)
  }
  if (sum(pair$win) <= 0 || sum(pair$loss) <= 0) {
    stop("Regression is unidentified because all pairwise time is one-sided.",
         call. = FALSE)
  }
  objective <- function(beta) {
    value <- .rmtif_regression_components(
      beta, pair$win, pair$loss, design$z0, design$z1
    )
    0.5 * sum(value$score^2)
  }
  gradient <- function(beta) {
    value <- .rmtif_regression_components(
      beta, pair$win, pair$loss, design$z0, design$z1
    )
    as.numeric(t(value$jacobian) %*% value$score)
  }
  fit <- stats::optim(
    par = as.numeric(beta_initial), fn = objective, gr = gradient,
    method = "BFGS",
    control = list(maxit = as.integer(maximum_iterations), reltol = 1e-11)
  )
  beta <- fit$par
  names(beta) <- design$names
  components <- .rmtif_regression_components(
    beta, pair$win, pair$loss, design$z0, design$z1
  )
  sensitivity_condition_number <- .rmtif_condition_number(components$jacobian)
  if (!is.finite(sensitivity_condition_number)) {
    stop("Regression score sensitivity is singular or nearly singular; ",
         "remove redundant covariates before reporting coefficients or intervals.",
         call. = FALSE)
  }
  n0 <- pair$split$n0
  n1 <- pair$split$n1
  phi0 <- matrix(0, n0, q)
  phi1 <- matrix(0, n1, q)
  row_residual <- rowSums(components$residual)
  column_residual <- colSums(components$residual)
  phi0 <- (components$residual %*% design$z1 -
             design$z0 * row_residual) / n1
  phi1 <- (design$z1 * column_residual -
             t(components$residual) %*% design$z0) / n0
  score_variance <- stats::cov(phi0) / n0 + stats::cov(phi1) / n1
  inverse_jacobian <- .rmtif_pinv(components$jacobian)
  vcov <- inverse_jacobian %*% score_variance %*% t(inverse_jacobian)
  rownames(vcov) <- colnames(vcov) <- design$names
  se <- sqrt(pmax(0, diag(vcov)))
  critical <- stats::qnorm(1 - (1 - conf_level) / 2)
  conf_int <- cbind(lower = beta - critical * se,
                    upper = beta + critical * se)
  statistic <- rep(NA_real_, length(beta))
  names(statistic) <- names(beta)
  has_positive_se <- is.finite(se) & se > 0
  statistic[has_positive_se] <- beta[has_positive_se] / se[has_positive_se]
  p_value <- rep(NA_real_, length(beta))
  names(p_value) <- names(beta)
  has_statistic <- is.finite(statistic)
  p_value[has_statistic] <- 2 * stats::pnorm(-abs(statistic[has_statistic]))
  score_norm <- sqrt(sum(components$score^2))
  root_converged <- isTRUE(fit$convergence == 0L) &&
    is.finite(score_norm) && score_norm <= score_tolerance &&
    !components$exponent_clipped
  list(
    coefficients = beta,
    vcov = vcov,
    se = se,
    conf_int = conf_int,
    statistic = statistic,
    p_value = p_value,
    score = components$score,
    score_norm = score_norm,
    score_tolerance = score_tolerance,
    convergence = fit$convergence,
    root_converged = root_converged,
    message = fit$message,
    objective = fit$value,
    exponent_clipped = components$exponent_clipped,
    maximum_abs_exponent = components$maximum_abs_exponent,
    sensitivity_condition_number = sensitivity_condition_number,
    inference = "two-sample projection"
  )
}

#' Fit covariate-adjusted pairwise rmt-IF regression
#'
#' The model uses `treatment + (X_treatment - X_control)` as its pairwise
#' design, so the first coefficient is the treatment effect for equal baseline
#' covariates.  The score is exactly the manuscript score
#' `(win - loss * exp(beta^T difference)) * difference`.  With arm-specific
#' IPCW, a stratified bootstrap refits both censoring curves and is the primary
#' inference; the fixed-weight sandwich is retained as a diagnostic. A singular
#' or nearly singular score-sensitivity matrix makes coefficients and intervals
#' unidentified: the function raises an error instead of reporting a
#' pseudoinverse-based, spuriously narrow interval. Near singular means the
#' smallest singular value is at most the largest times square-root machine
#' precision; remove redundant covariates and refit.
#'
#' @param subjects List of `rmtif_subject` objects.
#' @param tau Restriction horizon.
#' @param censoring Either `"none"` or terminal-aware `"arm_km"`.
#' @param bootstrap Stratified bootstrap replicates. `0` returns fixed-weight
#'   diagnostic inference only.
#' @param conf_level Confidence level.
#' @param bootstrap_success_fraction Minimum fraction of finite bootstrap
#'   fits needed to report bootstrap inference.
#' @param seed Optional reproducibility seed.
#' @param minimum_survival Positivity threshold for each censoring curve.
#' @param beta_initial Optional starting coefficient vector.
#' @return An `rmtif_regression` fit with treatment and covariate coefficients.
#' @export
rmtif_regression <- function(subjects, tau,
                             censoring = c("none", "arm_km"),
                             bootstrap, conf_level,
                             bootstrap_success_fraction, seed = NULL,
                             minimum_survival = 0, beta_initial = NULL) {
  censoring <- match.arg(censoring)
  conf_level <- .rmtif_conf_level(conf_level)
  bootstrap_success_fraction <- .rmtif_bootstrap_success_fraction(
    bootstrap_success_fraction
  )
  minimum_survival <- .rmtif_minimum_survival(minimum_survival)
  bootstrap <- .rmtif_scalar_count(bootstrap, "bootstrap", lower = 0)
  if (bootstrap == 1L) {
    stop("bootstrap must be 0 or at least 2.", call. = FALSE)
  }
  pair <- .rmtif_pair_matrices(subjects, tau, censoring, minimum_survival)
  fixed <- .rmtif_regression_from_pair(pair, conf_level, beta_initial)
  if (!isTRUE(fixed$root_converged)) {
    stop(
      "The regression optimizer did not produce a verified score root ",
      "(convergence = ", fixed$convergence,
      ", score norm = ", format(fixed$score_norm, digits = 6),
      ", exponent clipped = ", fixed$exponent_clipped, ").",
      call. = FALSE
    )
  }
  fixed$fixed_weight_vcov <- fixed$vcov
  fixed$fixed_weight_se <- fixed$se
  fixed$bootstrap_coefficients <- matrix(numeric(), nrow = 0L,
                                         ncol = length(fixed$coefficients))
  fixed$bootstrap_success <- 0L
  fixed$censoring <- censoring
  fixed$n0 <- pair$split$n0
  fixed$n1 <- pair$split$n1
  if (bootstrap >= 2L) {
    split <- pair$split
    boot <- .rmtif_with_seed(seed, {
      result <- matrix(NA_real_, nrow = bootstrap,
                       ncol = length(fixed$coefficients))
      for (b in seq_len(bootstrap)) {
        sample <- c(
          split$control[sample.int(split$n0, split$n0, replace = TRUE)],
          split$treatment[sample.int(split$n1, split$n1, replace = TRUE)]
        )
        result[b, ] <- tryCatch({
          pair_b <- .rmtif_pair_matrices(sample, tau, censoring,
                                         minimum_survival)
          fitted_b <- .rmtif_regression_from_pair(
            pair_b, conf_level, beta_initial = fixed$coefficients
          )
          if (!isTRUE(fitted_b$root_converged)) {
            stop("Bootstrap regression fit did not produce a verified root.",
                 call. = FALSE)
          }
          fitted_b$coefficients
        }, error = function(e) rep(NA_real_, ncol(result)))
      }
      result
    })
    complete <- stats::complete.cases(boot)
    if (sum(complete) < max(2L, ceiling(
      bootstrap_success_fraction * bootstrap
    ))) {
      stop("Too few regression bootstrap fits succeeded.",
           call. = FALSE)
    }
    good <- boot[complete, , drop = FALSE]
    colnames(good) <- names(fixed$coefficients)
    fixed$bootstrap_coefficients <- boot
    fixed$bootstrap_success <- nrow(good)
    fixed$vcov <- stats::cov(good)
    fixed$se <- sqrt(diag(fixed$vcov))
    critical <- stats::qnorm(1 - (1 - conf_level) / 2)
    fixed$conf_int <- cbind(
      lower = fixed$coefficients - critical * fixed$se,
      upper = fixed$coefficients + critical * fixed$se
    )
    fixed$statistic <- rep(NA_real_, length(fixed$coefficients))
    names(fixed$statistic) <- names(fixed$coefficients)
    has_positive_se <- is.finite(fixed$se) & fixed$se > 0
    fixed$statistic[has_positive_se] <-
      fixed$coefficients[has_positive_se] / fixed$se[has_positive_se]
    fixed$p_value <- rep(NA_real_, length(fixed$coefficients))
    names(fixed$p_value) <- names(fixed$coefficients)
    has_statistic <- is.finite(fixed$statistic)
    fixed$p_value[has_statistic] <-
      2 * stats::pnorm(-abs(fixed$statistic[has_statistic]))
    fixed$inference <- paste0(
      "treatment-stratified bootstrap (", nrow(good), "/", bootstrap,
      " successful; censoring curves refit)"
    )
  } else if (censoring == "arm_km") {
    fixed$inference <- "fixed-weight projection diagnostic"
  }
  class(fixed) <- "rmtif_regression"
  fixed
}


# Poisson attempted-event updates for reversible ordered-state processes
# -----------------------------------------------------------------------------

#' Convert a Poisson attempted-event update model to a CTMC generator
#'
#' Attempted updates occur at the state-specific Poisson rates. Each attempt
#' draws a destination from the corresponding row of `transition`. A diagonal
#' draw leaves the state unchanged; it is not a state transition. The resulting
#' generator has off-diagonal entries `rates[i] * transition[i, j]` and rows
#' summing to zero. This construction accommodates reversible movement and
#' positive diagonal update probabilities.
#'
#' @param transition Row-stochastic attempted-update probability matrix.
#' @param rates Positive attempted-update rates, one per state.
#' @return CTMC generator matrix with the same row and column names.
#' @export
rmtif_poisson_update_generator <- function(transition, rates) {
  transition <- as.matrix(transition)
  rates <- as.numeric(rates)
  if (!is.numeric(transition) || nrow(transition) != ncol(transition) ||
      nrow(transition) < 2L || length(rates) != nrow(transition) ||
      anyNA(transition) || any(!is.finite(transition)) ||
      any(transition < 0) || any(abs(rowSums(transition) - 1) > 1e-8) ||
      anyNA(rates) || any(!is.finite(rates)) || any(rates <= 0)) {
    stop("transition must be row-stochastic and rates must be positive.",
         call. = FALSE)
  }
  generator <- sweep(transition, 1L, rates, `*`)
  diag(generator) <- 0
  diag(generator) <- -rowSums(generator)
  generator
}

#' Simulate finite-state Poisson attempted-event update trajectories
#'
#' At each state an exponential waiting time with the specified attempted-event
#' rate precedes a destination draw from `transition`. A self-update advances
#' time while retaining the same state. The output is one complete observed
#' ordered-state trajectory per subject; the supplied state order determines
#' favorable versus unfavorable states in later rmt-IF calculations.
#'
#' @param n Number of subjects.
#' @param transition Row-stochastic attempted-update probability matrix.
#' @param rates Positive attempted-update rates, one per state.
#' @param tau Follow-up horizon.
#' @param states Ordered numeric state labels.
#' @param initial_state Initial state label.
#' @param arm Treatment arm, 0 or 1.
#' @param absorbing_states Optional explicit absorbing state labels.
#' @param seed Optional seed.
#' @return List of complete `rmtif_subject` trajectories.
#' @export
rmtif_simulate_poisson_updates <- function(n, transition, rates, tau,
                                           states, initial_state, arm,
                                           absorbing_states, seed = NULL) {
  n <- .rmtif_scalar_count(n, "n", lower = 1)
  tau <- .rmtif_scalar_numeric(tau, "tau", lower = 0)
  transition <- as.matrix(transition)
  rates <- as.numeric(rates)
  states <- as.numeric(states)
  if (!is.numeric(transition) || !length(states) ||
      anyNA(states) || any(!is.finite(states)) || anyDuplicated(states) ||
      nrow(transition) != ncol(transition) ||
      nrow(transition) != length(states) ||
      length(rates) != length(states) || anyNA(rates) ||
      any(!is.finite(rates)) || any(rates <= 0) ||
      anyNA(transition) || any(!is.finite(transition)) ||
      any(transition < 0) || any(abs(rowSums(transition) - 1) > 1e-8)) {
    stop("transition, rates, and states have incompatible or invalid values.",
         call. = FALSE)
  }
  if (length(initial_state) != 1L || is.na(initial_state) ||
      !initial_state %in% states || !all(absorbing_states %in% states)) {
    stop("initial_state and absorbing_states must be in states.", call. = FALSE)
  }
  .rmtif_with_seed(seed, {
    result <- vector("list", n)
    for (id in seq_len(n)) {
      time <- 0
      state <- initial_state
      times <- 0
      path <- state
      while (time < tau && !state %in% absorbing_states) {
        index <- match(state, states)
        next_time <- time + stats::rexp(1L, rates[index])
        if (next_time >= tau) {
          break
        }
        state <- sample(states, 1L, prob = transition[index, ])
        time <- next_time
        times <- c(times, time)
        path <- c(path, state)
      }
      result[[id]] <- rmtif_subject(
        arm = arm, times = times, states = path, tau = tau,
        censor_time = tau, censor_event = FALSE, id = id
      )
    }
    result
  })
}

#' Calibrate complete-path projection inference under Poisson updates
#'
#' Repeatedly simulate two arms of a reversible ordered-state process and
#' compute the complete-path rmt-IF curve at each requested horizon. The
#' population truth follows from the corresponding CTMC generators. The
#' returned summary reports bias, empirical SD, projection SE, coverage and
#' a two-sided zero-contrast rejection proportion for every horizon.
#'
#' @param scenario Either `"Null"` or an alternative design.
#' @param n_per_arm Subjects per arm.
#' @param replications Monte Carlo replications.
#' @param tau Restricted horizon.
#' @param increment Evaluation-grid spacing.
#' @param seed Reproducibility seed.
#' @param design List with `transition_control`, `transition_treatment`,
#'   `rates_control`, `rates_treatment`, `states`, `initial_state`, and
#'   `absorbing_states`; all study values are supplied by the caller.
#' @param conf_level Confidence level.
#' @param test_level Rejection threshold.
#' @param truth_grid_points Number of numerical integration points for the
#'   population CTMC contrast.
#' @return List with replicate-level curves and a pointwise summary containing
#'   population truth, bias, empirical SD, mean projection SE, normal-interval
#'   coverage, and zero-null rejection proportion.
#' @export
rmtif_poisson_update_study <- function(
    scenario, n_per_arm, replications, tau, increment, seed, design,
    conf_level, test_level, truth_grid_points) {
  n_per_arm <- .rmtif_scalar_count(
    n_per_arm, "n_per_arm", lower = 2
  )
  replications <- .rmtif_scalar_count(
    replications, "replications", lower = 1
  )
  increment <- .rmtif_scalar_numeric(increment, "increment", lower = 0)
  if (increment <= 0 || increment > tau) {
    stop("increment must be in (0, tau].", call. = FALSE)
  }
  transition_control <- design$transition_control
  transition_treatment <- if (scenario == "Null") {
    transition_control
  } else {
    design$transition_treatment
  }
  rates_control <- design$rates_control
  rates_treatment <- if (scenario == "Null") {
    rates_control
  } else {
    design$rates_treatment
  }
  generator_control <- rmtif_poisson_update_generator(
    transition_control, rates_control
  )
  generator_treatment <- rmtif_poisson_update_generator(
    transition_treatment, rates_treatment
  )
  evaluation_times <- seq(increment, tau, by = increment)
  if (utils::tail(evaluation_times, 1L) < tau) {
    evaluation_times <- c(evaluation_times, tau)
  }
  truth_by_time <- if (scenario == "Null") {
    rep(0, length(evaluation_times))
  } else {
    vapply(evaluation_times, function(horizon) {
      rmtif_true_ctmc_contrast(
        generator_control, generator_treatment,
        initial_state = design$initial_state, tau = horizon,
        grid_points = truth_grid_points
      )
    }, numeric(1))
  }
  .rmtif_with_seed(seed, {
    seeds <- sample.int(.Machine$integer.max, 2L * replications)
    curves <- vector("list", replications)
    for (replicate in seq_len(replications)) {
      control <- rmtif_simulate_poisson_updates(
        n = n_per_arm, transition = transition_control,
        rates = rates_control, tau = tau, states = design$states,
        initial_state = design$initial_state, arm = 0,
        absorbing_states = design$absorbing_states,
        seed = seeds[2L * replicate - 1L]
      )
      treatment <- rmtif_simulate_poisson_updates(
        n = n_per_arm, transition = transition_treatment,
        rates = rates_treatment, tau = tau, states = design$states,
        initial_state = design$initial_state, arm = 1,
        absorbing_states = design$absorbing_states,
        seed = seeds[2L * replicate]
      )
      curve <- rmtif_complete_state_count_curve(
        c(control, treatment), evaluation_times, conf_level = conf_level
      )
      curve$true <- truth_by_time
      curve$replicate <- replicate
      curves[[replicate]] <- curve
    }
    replicates <- do.call(rbind, curves)
    rownames(replicates) <- NULL
    summary <- do.call(rbind, lapply(split(replicates, replicates$time),
                                     function(x) {
      finite_estimate <- is.finite(x$estimate)
      finite_inference <- finite_estimate & is.finite(x$se) &
        is.finite(x$lower) & is.finite(x$upper)
      finite_test <- finite_inference & is.finite(x$p_value)
      truth <- x$true[1L]
      data.frame(
        scenario = scenario,
        time = x$time[1L],
        true = truth,
        mean_estimate = mean(x$estimate[finite_estimate]),
        bias = mean(x$estimate[finite_estimate]) - truth,
        empirical_sd = stats::sd(x$estimate[finite_estimate]),
        mean_projection_se = mean(x$se[finite_inference]),
        coverage = mean(x$lower[finite_inference] <= truth &
                             x$upper[finite_inference] >= truth),
        rejection = mean(x$p_value[finite_test] < test_level),
        empirical_variance = stats::var(x$estimate[finite_estimate]),
        mean_projection_variance = mean(x$se[finite_inference]^2),
        monte_carlo_lower = mean(x$estimate[finite_estimate]) -
          stats::qnorm(1 - (1 - conf_level) / 2) *
          stats::sd(x$estimate[finite_estimate]) /
          sqrt(sum(finite_estimate)),
        monte_carlo_upper = mean(x$estimate[finite_estimate]) +
          stats::qnorm(1 - (1 - conf_level) / 2) *
          stats::sd(x$estimate[finite_estimate]) /
          sqrt(sum(finite_estimate)),
        successful = sum(finite_inference),
        replications = nrow(x), row.names = NULL
      )
    }))
    rownames(summary) <- NULL
    list(replicates = replicates, summary = summary)
  })
}

# -----------------------------------------------------------------------------
# Exactly specified Gumbel pairwise-regression sanity check
# -----------------------------------------------------------------------------

.rmtif_rgumbel <- function(n) {
  -log(-log(stats::runif(n)))
}

#' Simulate a Gumbel pairwise-regression data-generating mechanism
#'
#' Generate Gaussian baseline covariates with the supplied covariance matrix,
#' treatment allocation, and independent standard Gumbel errors. The resulting
#' continuous latent utility is constant over follow-up and is represented by
#' its negative as an ordered state, so smaller states mean better outcomes.
#' This is a regression-equation sanity check rather than a multistate DGM.
#'
#' @param n_per_arm Subjects per arm.
#' @param beta_true Named vector: treatment intercept followed by slopes.
#' @param tau Restriction horizon.
#' @param sigma_x Positive-definite covariance matrix for baseline covariates.
#' @param seed Optional seed.
#' @return Complete `rmtif_subject` list and true parameter.
#' @export
rmtif_simulate_gumbel_trial <- function(n_per_arm, beta_true, tau,
                                        sigma_x, seed = NULL) {
  n_per_arm <- .rmtif_scalar_count(
    n_per_arm, "n_per_arm", lower = 2
  )
  beta_true <- as.numeric(beta_true)
  names(beta_true) <- c("treatment", if (length(beta_true) > 1L) {
    paste0("covariate", seq_len(length(beta_true) - 1L))
  } else character())
  p <- length(beta_true) - 1L
  if (p < 1L) {
    stop("beta_true must contain a treatment effect and at least one slope.",
         call. = FALSE)
  }
  sigma_x <- as.matrix(sigma_x)
  if (!all(dim(sigma_x) == c(p, p))) {
    stop("sigma_x must be p by p.", call. = FALSE)
  }
  .rmtif_with_seed(seed, {
    n <- 2L * n_per_arm
    arm <- rep(0:1, each = n_per_arm)
    x <- matrix(stats::rnorm(n * p), nrow = n) %*% chol(sigma_x)
    utility <- beta_true[1L] * arm + as.numeric(x %*% beta_true[-1L]) +
      .rmtif_rgumbel(n)
    subjects <- lapply(seq_len(n), function(i) {
      rmtif_subject(
        arm = arm[i], times = 0, states = -utility[i], tau = tau,
        covariates = x[i, ], censor_time = tau, id = i
      )
    })
    list(subjects = subjects, beta_true = beta_true, tau = tau)
  })
}

#' Simulate the constant-outcome Gumbel regression DGM
#'
#' This descriptive entry point generates the same Gaussian-covariate,
#' Gumbel-error pairwise-regression trial as `rmtif_simulate_gumbel_trial`.
#' All treatment, slope, covariance, horizon, and sample-size inputs come from
#' the caller.
#'
#' @param n_per_arm Subjects per treatment arm.
#' @param beta_true Named treatment-intercept and covariate-slope vector.
#' @param tau Restricted time horizon.
#' @param sigma_x Positive-definite baseline-covariate covariance matrix.
#' @param seed Optional reproducibility seed.
#' @return List with complete ordered-state subjects, true parameter vector,
#'   and restriction horizon.
#' @export
rmtif_simulate_gumbel_regression <- function(n_per_arm, beta_true, tau,
                                              sigma_x, seed = NULL) {
  rmtif_simulate_gumbel_trial(n_per_arm, beta_true, tau, sigma_x, seed)
}

.rmtif_gumbel_fit_fast <- function(simulation, conf_level) {
  conf_level <- .rmtif_conf_level(conf_level)
  split <- .rmtif_validate_subjects(simulation$subjects, simulation$tau,
                                    minimum_per_arm = 2L)
  design <- .rmtif_regression_design(list(split = split))
  utility0 <- -vapply(split$control, function(x) x$states[1L], numeric(1))
  utility1 <- -vapply(split$treatment, function(x) x$states[1L], numeric(1))
  order0 <- order(utility0)
  utility0_sorted <- utility0[order0]
  z0 <- design$z0
  z0_sorted <- z0[order0, , drop = FALSE]
  z1 <- design$z1
  n0 <- split$n0
  n1 <- split$n1
  q <- ncol(z0)

  # The Gumbel process is constant in time. Sorting the controls reduces each
  # score/Jacobian evaluation from O(n0*n1) dense pair operations to cumulative
  # O((n0+n1)*q^2) operations. The final sandwich still uses the exact dense
  # pair residual once, so this is algebraically identical to the general fit.
  prefix_z <- apply(z0_sorted, 2L, cumsum)
  if (q == 1L) {
    prefix_z <- matrix(prefix_z, ncol = 1L)
  }
  prefix_z_padded <- rbind(rep(0, q), prefix_z)
  below_strict <- findInterval(utility1, utility0_sorted,
                               left.open = TRUE)
  first_above <- findInterval(utility1, utility0_sorted) + 1L
  win_sum <- colSums(z1 * below_strict -
                       prefix_z_padded[below_strict + 1L, , drop = FALSE])
  fast_components <- function(beta) {
    eta0 <- as.numeric(z0_sorted %*% beta)
    eta1 <- as.numeric(z1 %*% beta)
    weight0 <- exp(pmax(-700, pmin(700, -eta0)))
    suffix_weight <- c(rev(cumsum(rev(weight0))), 0)
    suffix_z <- apply(z0_sorted * weight0, 2L,
                      function(x) rev(cumsum(rev(x))))
    if (q == 1L) suffix_z <- matrix(suffix_z, ncol = 1L)
    suffix_z <- rbind(suffix_z, rep(0, q))
    suffix_outer <- array(0, dim = c(n0 + 1L, q, q))
    for (r in seq_len(q)) {
      for (s in seq_len(q)) {
        value <- weight0 * z0_sorted[, r] * z0_sorted[, s]
        suffix_outer[seq_len(n0), r, s] <- rev(cumsum(rev(value)))
      }
    }
    multiplier <- exp(pmax(-700, pmin(700, eta1)))
    selected_weight <- suffix_weight[first_above]
    selected_z <- suffix_z[first_above, , drop = FALSE]
    loss_sum <- colSums((z1 * selected_weight - selected_z) * multiplier)
    term11 <- crossprod(z1, z1 * (multiplier * selected_weight))
    term10 <- crossprod(z1, selected_z * multiplier)
    term00 <- matrix(0, q, q)
    for (r in seq_len(q)) {
      for (s in seq_len(q)) {
        term00[r, s] <- sum(multiplier * suffix_outer[first_above, r, s])
      }
    }
    scale <- simulation$tau / (n0 * n1)
    list(score = as.numeric((win_sum - loss_sum) * scale),
         jacobian = -(term11 - term10 - t(term10) + term00) * scale)
  }
  objective <- function(beta) {
    value <- fast_components(beta)$score
    0.5 * sum(value^2)
  }
  gradient <- function(beta) {
    value <- fast_components(beta)
    as.numeric(t(value$jacobian) %*% value$score)
  }
  optimized <- stats::optim(
    numeric(q), objective, gradient, method = "BFGS",
    control = list(maxit = 300L, reltol = 1e-11)
  )

  state0 <- -utility0
  state1 <- -utility1
  pair <- list(
    win = simulation$tau * (outer(state0, state1, ">") + 0),
    loss = simulation$tau * (outer(state0, state1, "<") + 0),
    split = split
  )
  components <- .rmtif_regression_components(
    optimized$par, pair$win, pair$loss, design$z0, design$z1
  )
  sensitivity_condition_number <- .rmtif_condition_number(components$jacobian)
  if (!is.finite(sensitivity_condition_number)) {
    stop("Regression score sensitivity is singular or nearly singular; ",
         "remove redundant covariates before reporting coefficients or intervals.",
         call. = FALSE)
  }
  row_residual <- rowSums(components$residual)
  column_residual <- colSums(components$residual)
  phi0 <- (components$residual %*% design$z1 -
             design$z0 * row_residual) / n1
  phi1 <- (design$z1 * column_residual -
             t(components$residual) %*% design$z0) / n0
  score_variance <- stats::cov(phi0) / n0 + stats::cov(phi1) / n1
  inverse_jacobian <- .rmtif_pinv(components$jacobian)
  vcov <- inverse_jacobian %*% score_variance %*% t(inverse_jacobian)
  rownames(vcov) <- colnames(vcov) <- design$names
  coefficients <- optimized$par
  names(coefficients) <- design$names
  se <- sqrt(pmax(0, diag(vcov)))
  critical <- stats::qnorm(1 - (1 - conf_level) / 2)
  conf_int <- cbind(
    lower = coefficients - critical * se,
    upper = coefficients + critical * se
  )
  statistic <- coefficients / se
  score_norm <- sqrt(sum(components$score^2))
  score_tolerance <- 1e-7
  root_converged <- isTRUE(optimized$convergence == 0L) &&
    is.finite(score_norm) && score_norm <= score_tolerance &&
    !components$exponent_clipped
  list(
    coefficients = coefficients, vcov = vcov, se = se,
    conf_int = conf_int, statistic = statistic,
    p_value = 2 * stats::pnorm(-abs(statistic)),
    score = components$score, score_norm = score_norm,
    score_tolerance = score_tolerance,
    convergence = optimized$convergence, root_converged = root_converged,
    exponent_clipped = components$exponent_clipped,
    maximum_abs_exponent = components$maximum_abs_exponent,
    message = optimized$message, objective = optimized$value,
    inference = "two-sample projection"
  )
}

#' Summarize Gumbel regression recovery over simulated trials
#'
#' The summary contains the mean estimate, empirical SD, mean estimated SE,
#' interval coverage, and convergence counts for every parameter and sample
#' size. Trials with singular or nearly singular sensitivity are excluded from
#' summary estimates and counted as failed fits.
#'
#' @param n_per_arm_values Per-arm sample sizes.
#' @param replications Monte Carlo trials per sample size.
#' @param beta_true True treatment and covariate coefficients.
#' @param tau Restriction horizon.
#' @param sigma_x Baseline covariate covariance matrix.
#' @param conf_level Confidence level.
#' @param seed Reproducibility seed.
#' @return List containing replicate-level and summary data frames.
#' @export
rmtif_gumbel_recovery_study <- function(
    n_per_arm_values, replications, beta_true, tau, sigma_x,
    conf_level, seed) {
  replications <- .rmtif_scalar_count(
    replications, "replications", lower = 1
  )
  .rmtif_with_seed(seed, {
    seeds <- sample.int(.Machine$integer.max, length(n_per_arm_values) * replications)
    records <- vector("list", length(seeds) * length(beta_true))
    index <- 0L
    seed_index <- 0L
    for (n_per_arm in n_per_arm_values) {
      for (replicate in seq_len(replications)) {
        seed_index <- seed_index + 1L
        simulation <- rmtif_simulate_gumbel_regression(
          n_per_arm, beta_true, tau, sigma_x = sigma_x,
          seed = seeds[seed_index]
        )
        fit <- tryCatch(.rmtif_gumbel_fit_fast(simulation, conf_level),
                        error = function(e) NULL)
        for (parameter in seq_along(beta_true)) {
          index <- index + 1L
          records[[index]] <- data.frame(
            n_per_arm = n_per_arm,
            replicate = replicate,
            parameter = names(simulation$beta_true)[parameter],
            true = simulation$beta_true[parameter],
            estimate = if (is.null(fit)) NA_real_ else fit$coefficients[parameter],
            estimated_se = if (is.null(fit)) NA_real_ else fit$se[parameter],
            lower = if (is.null(fit)) NA_real_ else fit$conf_int[parameter, "lower"],
            upper = if (is.null(fit)) NA_real_ else fit$conf_int[parameter, "upper"],
            converged = !is.null(fit) && isTRUE(fit$root_converged)
          )
        }
      }
    }
    replicates <- do.call(rbind, records)
    groups <- split(replicates, list(replicates$n_per_arm,
                                     replicates$parameter), drop = TRUE)
    summary <- do.call(rbind, lapply(groups, function(x) {
      good <- x$converged & is.finite(x$estimate) & is.finite(x$estimated_se)
      data.frame(
        n_per_arm = x$n_per_arm[1L],
        parameter = x$parameter[1L],
        true = x$true[1L],
        mean_estimate = mean(x$estimate[good]),
        empirical_sd = stats::sd(x$estimate[good]),
        mean_estimated_se = mean(x$estimated_se[good]),
        coverage = mean(x$lower[good] <= x$true[good] &
                             x$upper[good] >= x$true[good]),
        successful = sum(good),
        replications = nrow(x),
        row.names = NULL
      )
    }))
    rownames(summary) <- NULL
    summary <- summary[
      order(summary$n_per_arm, match(summary$parameter, names(beta_true))),
      , drop = FALSE
    ]
    rownames(summary) <- NULL
    list(replicates = replicates, summary = summary)
  })
}

# -----------------------------------------------------------------------------
# Continuous-time Markov machinery and dedicated differential-censoring study
# -----------------------------------------------------------------------------

#' Build a continuous-time Markov generator from transition intensities
#'
#' The caller supplies nonnegative off-diagonal transition intensities in a
#' square matrix. Diagonal values are ignored and replaced by the negative
#' off-diagonal row sums. This permits reversible movement among nonterminal
#' states and an absorbing terminal state when its off-diagonal row is zero.
#'
#' @param off_diagonal_rates Square matrix of transition intensities.
#' @return Valid CTMC generator with zero row sums.
#' @export
rmtif_ctmc_generator <- function(off_diagonal_rates) {
  rates <- as.matrix(off_diagonal_rates)
  if (!is.numeric(rates) || nrow(rates) != ncol(rates) ||
      nrow(rates) < 2L || anyNA(rates) || any(!is.finite(rates))) {
    stop("off_diagonal_rates must be a finite square numeric matrix.",
         call. = FALSE)
  }
  diag(rates) <- 0
  if (any(rates < 0)) {
    stop("Off-diagonal transition intensities must be nonnegative.",
         call. = FALSE)
  }
  diag(rates) <- -rowSums(rates)
  .rmtif_validate_generator(rates)
}

.rmtif_validate_generator <- function(generator) {
  generator <- as.matrix(generator)
  if (!is.numeric(generator) || nrow(generator) != ncol(generator) ||
      nrow(generator) < 2L || any(!is.finite(generator))) {
    stop("generator must be a finite square numeric matrix.", call. = FALSE)
  }
  diagonal <- diag(generator)
  off_diagonal <- generator
  diag(off_diagonal) <- 0
  if (any(off_diagonal < 0) || any(diagonal > 1e-12) ||
      any(abs(rowSums(generator)) > 1e-8)) {
    stop("generator must have nonnegative off-diagonals and zero row sums.",
         call. = FALSE)
  }
  generator
}

#' Simulate a complete continuous-time Markov state path
#'
#' Wait at each state according to its generator exit intensity, then choose
#' the next state in proportion to the off-diagonal intensities. An absorbing
#' state ends transitions while retaining its state through the horizon.
#'
#' @param generator Valid CTMC generator matrix.
#' @param initial_state Initial state index starting at zero.
#' @param tau Follow-up horizon.
#' @return List with right-continuous `times` and `states` vectors.
#' @export
rmtif_simulate_ctmc_path <- function(generator, initial_state, tau) {
  generator <- .rmtif_validate_generator(generator)
  initial_state <- .rmtif_scalar_count(initial_state, "initial_state")
  tau <- .rmtif_scalar_numeric(tau, "tau", lower = 0)
  state <- initial_state + 1L
  if (state < 1L || state > nrow(generator)) {
    stop("initial_state is outside the generator state space.", call. = FALSE)
  }
  time <- 0
  times <- 0
  states <- state - 1L
  while (time < tau) {
    rate <- -generator[state, state]
    if (rate <= 0) {
      break
    }
    next_time <- time + stats::rexp(1L, rate)
    if (next_time >= tau) {
      break
    }
    probabilities <- generator[state, ]
    probabilities[state] <- 0
    probabilities <- probabilities / sum(probabilities)
    state <- sample.int(nrow(generator), 1L, prob = probabilities)
    time <- next_time
    times <- c(times, time)
    states <- c(states, state - 1L)
  }
  list(times = times, states = states)
}

.rmtif_terminal_time <- function(path, death_state) {
  index <- which(path$states == death_state)
  if (length(index)) path$times[index[1L]] else Inf
}

#' Apply loss to follow-up and terminal-event observation to a CTMC path
#'
#' A terminal state reached before censoring is observed exactly and carried
#' forward. Earlier loss to follow-up truncates the observed nonterminal path.
#'
#' @param path Complete CTMC path from `rmtif_simulate_ctmc_path`.
#' @param arm Treatment indicator.
#' @param tau Restriction horizon.
#' @param censor_time Independent loss-to-follow-up time.
#' @param death_state Absorbing terminal-state label.
#' @param covariates Optional numeric baseline covariate vector.
#' @param id Optional subject identifier.
#' @return One observed `rmtif_subject`.
#' @export
rmtif_observe_ctmc_path <- function(path, arm, tau, censor_time,
                                    death_state, covariates = numeric(),
                                    id = NULL) {
  latent_death <- .rmtif_terminal_time(path, death_state)
  observed_death <- if (latent_death <= min(censor_time, tau)) {
    latent_death
  } else {
    Inf
  }
  loss_event <- censor_time < min(latent_death, tau)
  observation_limit <- min(censor_time, observed_death, tau)
  keep <- path$times <= observation_limit
  rmtif_subject(
    arm = arm,
    times = path$times[keep],
    states = path$states[keep],
    tau = tau,
    covariates = covariates,
    censor_time = min(censor_time, tau),
    censor_event = loss_event,
    death_time = observed_death,
    death_state = if (is.finite(observed_death)) death_state else NULL,
    id = id
  )
}

#' Prepare control and treatment CTMC generators for a censoring study
#'
#' The design supplies every transition intensity and treatment multiplier.
#' Additional trial design fields are passed through for convenience.
#'
#' @param design List containing `rates_control`, `benefit_multiplier`,
#'   `initial_state`, `death_state`, `tau`, `primary_censoring_rates`,
#'   `allocations`, `trials`, and `bootstrap`.
#' @return Named list of control, null-treatment and benefit-treatment CTMC
#'   generators, state labels, and supplied trial settings.
#' @export
rmtif_ipcw_simulation_parameters <- function(design) {
  required <- c("rates_control", "benefit_multiplier", "initial_state",
                "death_state", "tau", "primary_censoring_rates",
                "allocations", "trials", "bootstrap")
  if (!is.list(design) || !all(required %in% names(design))) {
    stop("design must include: ", paste(required, collapse = ", "), ".",
         call. = FALSE)
  }
  rates <- design$rates_control
  multiplier <- design$benefit_multiplier
  if (!is.matrix(rates) || !is.matrix(multiplier) ||
      !identical(dim(rates), dim(multiplier)) ||
      !is.numeric(multiplier) || anyNA(multiplier) ||
      any(!is.finite(multiplier)) || any(multiplier < 0)) {
    stop("rates_control and benefit_multiplier must be equal-size finite ",
         "numeric matrices with nonnegative multipliers.", call. = FALSE)
  }
  list(
    generator_control = rmtif_ctmc_generator(rates),
    generator_treatment_null = rmtif_ctmc_generator(rates),
    generator_treatment_benefit = rmtif_ctmc_generator(rates * multiplier),
    initial_state = design$initial_state,
    death_state = design$death_state,
    tau = design$tau,
    censoring_rates = design$primary_censoring_rates,
    allocations = design$allocations,
    trials = design$trials,
    bootstrap = design$bootstrap
  )
}

#' Simulate one dedicated IPCW trial
#'
#' The aim is to evaluate the marginal arm-specific IPCW estimator under a
#' reversible nonterminal process, absorbing death, unequal allocation, and
#' differential but within-arm independent loss to follow-up.
#'
#' @param n0 Control sample size.
#' @param n1 Treatment sample size.
#' @param scenario `"Null"` or `"Treatment benefit"`.
#' @param tau Restriction horizon.
#' @param censoring_rates Arm-specific exponential loss-to-follow-up rates.
#' @param parameters Markov generators and state labels.
#' @param seed Optional seed.
#' @return List of observed subjects and the two population generators.
#' @export
rmtif_simulate_ipcw_trial <- function(
    n0, n1, scenario, tau, censoring_rates, parameters, seed = NULL) {
  n0 <- .rmtif_scalar_count(n0, "n0", lower = 2)
  n1 <- .rmtif_scalar_count(n1, "n1", lower = 2)
  if (!is.character(scenario) || length(scenario) != 1L ||
      !scenario %in% c("Null", "Treatment benefit")) {
    stop("scenario must be 'Null' or 'Treatment benefit'.", call. = FALSE)
  }
  generator0 <- parameters$generator_control
  generator1 <- if (scenario == "Null") {
    parameters$generator_treatment_null
  } else {
    parameters$generator_treatment_benefit
  }
  if (length(censoring_rates) != 2L || any(!is.finite(censoring_rates)) ||
      any(censoring_rates <= 0)) {
    stop("censoring_rates must contain two positive finite arm-specific rates.",
         call. = FALSE)
  }
  .rmtif_with_seed(seed, {
    arms <- c(rep(0L, n0), rep(1L, n1))
    loss_times <- c(stats::rexp(n0, censoring_rates[1L]),
                    stats::rexp(n1, censoring_rates[2L]))
    subjects <- lapply(seq_along(arms), function(i) {
      generator <- if (arms[i] == 0L) generator0 else generator1
      path <- rmtif_simulate_ctmc_path(generator, parameters$initial_state,
                                        tau)
      rmtif_observe_ctmc_path(
        path, arm = arms[i], tau = tau, censor_time = loss_times[i],
        death_state = parameters$death_state, id = i
      )
    })
    list(subjects = subjects, generator_control = generator0,
         generator_treatment = generator1, tau = tau, scenario = scenario)
  })
}

.rmtif_transition_probability <- function(generator, time,
                                          initial_state = 0L,
                                          tolerance = 1e-13,
                                          maximum_terms = 10000L) {
  generator <- .rmtif_validate_generator(generator)
  time <- .rmtif_scalar_numeric(time, "time", lower = 0)
  q <- max(-diag(generator))
  initial <- numeric(nrow(generator))
  initial[as.integer(initial_state) + 1L] <- 1
  if (time == 0 || q == 0) {
    return(initial)
  }
  qt <- q * time
  if (!is.finite(qt) || qt / 500 > maximum_terms) {
    stop("rate multiplied by time exceeds the supported uniformization range.",
         call. = FALSE)
  }
  if (qt > 500) {
    # Direct Poisson weights start at exp(-q*t), which underflows for large
    # horizons. Compute a well-scaled smaller-time transition matrix and
    # compose it using the Markov semigroup property.
    segments <- ceiling(qt / 500)
    step <- t(vapply(seq_len(nrow(generator)), function(state) {
      .rmtif_transition_probability(
        generator, time / segments, initial_state = state - 1L,
        tolerance = max(tolerance / segments, .Machine$double.eps),
        maximum_terms = maximum_terms
      )
    }, numeric(nrow(generator))))
    for (segment in seq_len(segments)) {
      initial <- as.numeric(initial %*% step)
    }
    return(pmax(0, initial / sum(initial)))
  }
  transition <- diag(nrow(generator)) + generator / q
  probability <- exp(-qt)
  state_probability <- initial
  result <- probability * state_probability
  cumulative <- probability
  k <- 0L
  while ((1 - cumulative) > tolerance && k < maximum_terms) {
    k <- k + 1L
    state_probability <- as.numeric(state_probability %*% transition)
    probability <- probability * qt / k
    result <- result + probability * state_probability
    cumulative <- cumulative + probability
  }
  if (k >= maximum_terms) {
    stop("Uniformization did not converge.", call. = FALSE)
  }
  pmax(0, result / sum(result))
}

#' Compute a population CTMC rmt-IF contrast
#'
#' The aim is to obtain simulation truth without Monte Carlo error.  Transition
#' probabilities are computed by uniformization and the pairwise state contrast
#' is integrated numerically.
#'
#' @param generator_control Control CTMC generator.
#' @param generator_treatment Treatment CTMC generator.
#' @param initial_state Common initial state, numbered from zero.
#' @param tau Restriction horizon.
#' @param grid_points Number of equally spaced integration points, supplied
#'   by the caller to control numerical accuracy.
#' @return Population treatment-minus-control rmt-IF in time units.
#' @export
rmtif_true_ctmc_contrast <- function(generator_control, generator_treatment,
                                     initial_state, tau,
                                     grid_points) {
  generator_control <- .rmtif_validate_generator(generator_control)
  generator_treatment <- .rmtif_validate_generator(generator_treatment)
  if (!all(dim(generator_control) == dim(generator_treatment))) {
    stop("The two generators must have identical dimensions.", call. = FALSE)
  }
  initial_state <- .rmtif_scalar_count(initial_state, "initial_state")
  if (initial_state >= nrow(generator_control)) {
    stop("initial_state is outside the generator state space.", call. = FALSE)
  }
  tau <- .rmtif_scalar_numeric(tau, "tau", lower = 0)
  grid_points <- .rmtif_scalar_count(grid_points, "grid_points", lower = 3)
  times <- seq(0, tau, length.out = grid_points)
  comparison <- outer(seq_len(nrow(generator_control)),
                      seq_len(nrow(generator_treatment)),
                      function(i, j) (j < i) - (j > i))
  value <- vapply(times, function(tt) {
    p0 <- .rmtif_transition_probability(generator_control, tt, initial_state)
    p1 <- .rmtif_transition_probability(generator_treatment, tt, initial_state)
    sum(outer(p0, p1) * comparison)
  }, numeric(1))
  width <- diff(times)
  sum(width * (utils::head(value, -1L) + utils::tail(value, -1L)) / 2)
}

#' Run the dedicated differential-censoring IPCW study
#'
#' The returned manuscript-ready summary contains truth, mean, bias, empirical
#' SD, mean bootstrap SE, normal-interval coverage, rejection probability,
#' and fit counts for every requested scenario/allocation combination.
#'
#' @param allocations List of named `c(n0, n1)` vectors.
#' @param scenarios Null and/or treatment-benefit scenarios.
#' @param trials Number of independently simulated trials.
#' @param bootstrap Stratified bootstrap replicates per trial.
#' @param tau Restriction horizon.
#' @param censoring_rates Named or ordered arm-specific exponential
#'   loss-to-follow-up rates used in every generated trial.
#' @param seed Reproducibility seed.
#' @param minimum_survival Positivity threshold.
#' @param parameters Markov generators and state labels.
#' @param conf_level Confidence level.
#' @param bootstrap_success_fraction Required successful bootstrap fraction.
#' @param test_level Rejection threshold.
#' @param truth_grid_points Number of numerical integration points for the
#'   population CTMC contrast.
#' @return List containing replicate-level and summary data frames.
#' @export
rmtif_ipcw_simulation_study <- function(
    allocations, scenarios, trials, bootstrap, tau, seed, minimum_survival,
    censoring_rates, parameters, conf_level, bootstrap_success_fraction,
    test_level, truth_grid_points) {
  trials <- .rmtif_scalar_count(trials, "trials", lower = 1)
  tau <- .rmtif_scalar_numeric(tau, "tau", lower = 0)
  bootstrap <- .rmtif_scalar_count(bootstrap, "bootstrap", lower = 2)
  minimum_survival <- .rmtif_minimum_survival(minimum_survival)
  if (length(censoring_rates) != 2L || any(!is.finite(censoring_rates)) ||
      any(censoring_rates <= 0)) {
    stop("censoring_rates must contain two positive finite arm-specific rates.",
         call. = FALSE)
  }
  allowed_scenarios <- c("Null", "Treatment benefit")
  if (!length(scenarios) || !all(scenarios %in% allowed_scenarios)) {
    stop("scenarios must contain 'Null' and/or 'Treatment benefit'.",
         call. = FALSE)
  }
  truth_by_scenario <- vapply(scenarios, function(scenario) {
    if (scenario == "Null") return(0)
    rmtif_true_ctmc_contrast(
      parameters$generator_control,
      parameters$generator_treatment_benefit,
      parameters$initial_state, tau, grid_points = truth_grid_points
    )
  }, numeric(1))
  names(truth_by_scenario) <- scenarios
  combinations <- expand.grid(
    allocation = seq_along(allocations), scenario = scenarios,
    stringsAsFactors = FALSE
  )
  .rmtif_with_seed(seed, {
    seeds <- sample.int(.Machine$integer.max, nrow(combinations) * trials)
    records <- vector("list", length(seeds))
    index <- 0L
    for (row in seq_len(nrow(combinations))) {
      allocation <- allocations[[combinations$allocation[row]]]
      scenario <- combinations$scenario[row]
      truth <- truth_by_scenario[[scenario]]
      for (replicate in seq_len(trials)) {
        index <- index + 1L
        trial <- rmtif_simulate_ipcw_trial(
          n0 = allocation["n0"], n1 = allocation["n1"],
          scenario = scenario, tau = tau,
          censoring_rates = censoring_rates, parameters = parameters,
          seed = seeds[index]
        )
        fit <- tryCatch(rmtif_ipcw(
          trial$subjects, tau, bootstrap = bootstrap,
          conf_level = conf_level,
          bootstrap_success_fraction = bootstrap_success_fraction,
          seed = .rmtif_next_seed(seeds[index]),
          minimum_survival = minimum_survival,
          projection_diagnostic = FALSE
        ), error = function(e) NULL)
        trial_arms <- vapply(trial$subjects, `[[`, integer(1), "arm")
        trial_ltfu <- vapply(trial$subjects, `[[`, logical(1), "censor_event")
        records[[index]] <- data.frame(
          scenario = scenario,
          n0 = as.integer(allocation["n0"]),
          n1 = as.integer(allocation["n1"]),
          replicate = replicate,
          true = truth,
          estimate = if (is.null(fit)) NA_real_ else fit$estimate,
          estimated_se = if (is.null(fit)) NA_real_ else fit$se,
          lower = if (is.null(fit)) NA_real_ else fit$conf_int["lower"],
          upper = if (is.null(fit)) NA_real_ else fit$conf_int["upper"],
          p_value = if (is.null(fit)) NA_real_ else fit$p_value,
          bootstrap_success = if (is.null(fit)) 0L else fit$bootstrap_success,
          ltfu_control = mean(trial_ltfu[trial_arms == 0L]),
          ltfu_treatment = mean(trial_ltfu[trial_arms == 1L])
        )
      }
    }
    replicates <- do.call(rbind, records)
    groups <- split(replicates,
                    list(replicates$scenario, replicates$n0, replicates$n1),
                    drop = TRUE)
    summary <- do.call(rbind, lapply(groups, function(x) {
      good <- is.finite(x$estimate) & is.finite(x$estimated_se)
      data.frame(
        scenario = x$scenario[1L], n0 = x$n0[1L], n1 = x$n1[1L],
        true = x$true[1L],
        mean_estimate = if (any(good)) mean(x$estimate[good]) else NA_real_,
        bias = if (any(good)) mean(x$estimate[good]) - x$true[1L] else NA_real_,
        empirical_sd = if (sum(good) >= 2L) stats::sd(x$estimate[good]) else NA_real_,
        mean_bootstrap_se = if (any(good)) mean(x$estimated_se[good]) else NA_real_,
        coverage = if (any(good)) {
          mean(x$lower[good] <= x$true[1L] &
                 x$upper[good] >= x$true[1L])
        } else NA_real_,
        rejection = if (any(good)) mean(x$p_value[good] < test_level) else NA_real_,
        mean_ltfu_control = mean(x$ltfu_control),
        mean_ltfu_treatment = mean(x$ltfu_treatment),
        successful = sum(good), trials = nrow(x),
        row.names = NULL
      )
    }))
    rownames(summary) <- NULL
    list(replicates = replicates, summary = summary)
  })
}

# -----------------------------------------------------------------------------
# Public Wcompo::hfmock illustration
# -----------------------------------------------------------------------------

#' Convert `Wcompo::hfmock` event records to scheduled-visit state paths
#'
#' At caller-specified visits, the better nonterminal state means no
#' hospitalization since the prior visit, the worse nonterminal state means at
#' least one hospitalization, and death is an exact absorbing state. An
#' end-of-record observation between visits censors the scheduled-visit process
#' at the most recent completed visit. The raw data do not distinguish every
#' such observation as loss to follow-up versus administrative closure.
#' A record extending to or beyond the requested horizon supports carrying
#' the most recent scheduled state to that horizon, including a horizon
#' between scheduled visits.
#'
#' @param data A data frame with `id`, `time`, `status`, `Training`, and
#'   `HF.etiology`, such as `Wcompo::hfmock`.
#' @param tau Restriction horizon.
#' @param visit_interval Scheduled-visit spacing.
#' @param raw_status Named vector with `censor`, `death`, and
#'   `hospitalization` raw event codes.
#' @param ordered_state Named vector with `no_hospitalization`,
#'   `hospitalization`, and `death` ordered-state codes.
#' @param quantile_probs Three requested probabilities for the arm audit's
#'   lower, median, and upper censor-time summaries.
#' @return List with `rmtif_subject` objects and raw-data audit summaries.
#' @export
rmtif_hfmock_subjects <- function(data, tau, visit_interval,
                                  raw_status, ordered_state, quantile_probs) {
  tau <- .rmtif_scalar_numeric(tau, "tau", lower = 0)
  visit_interval <- .rmtif_scalar_numeric(
    visit_interval, "visit_interval", lower = 0
  )
  if (tau <= 0 || visit_interval <= 0) {
    stop("tau and visit_interval must both be strictly positive.",
         call. = FALSE)
  }
  if (!is.numeric(quantile_probs) || length(quantile_probs) != 3L ||
      anyNA(quantile_probs) || any(quantile_probs < 0 | quantile_probs > 1) ||
      is.unsorted(quantile_probs, strictly = TRUE)) {
    stop("quantile_probs must be three increasing probabilities in [0, 1].",
         call. = FALSE)
  }
  if (!is.numeric(raw_status) ||
      !all(c("censor", "death", "hospitalization") %in% names(raw_status)) ||
      anyNA(raw_status) || anyDuplicated(unname(raw_status))) {
    stop("raw_status must contain distinct named numeric event codes.",
         call. = FALSE)
  }
  if (!is.numeric(ordered_state) ||
      !all(c("no_hospitalization", "hospitalization", "death") %in%
           names(ordered_state)) ||
      anyNA(ordered_state) ||
      !(ordered_state[["no_hospitalization"]] <
        ordered_state[["hospitalization"]] &&
        ordered_state[["hospitalization"]] < ordered_state[["death"]])) {
    stop("ordered_state must provide strictly increasing named state codes.",
         call. = FALSE)
  }
  required <- c("id", "time", "status", "Training", "HF.etiology")
  if (!is.data.frame(data) || !all(required %in% names(data))) {
    stop("data must contain: ", paste(required, collapse = ", "), ".",
         call. = FALSE)
  }
  if (anyNA(data[, required])) {
    stop("Required hfmock columns may not contain missing values.",
         call. = FALSE)
  }
  if (!is.numeric(data$time) || any(!is.finite(data$time)) ||
      any(data$time < 0)) {
    stop("hfmock event times must be finite nonnegative numbers.",
         call. = FALSE)
  }
  ids <- unique(data$id)
  subjects <- vector("list", length(ids))
  subject_summary <- vector("list", length(ids))
  for (i in seq_along(ids)) {
    records <- data[data$id == ids[i], , drop = FALSE]
    arm <- unique(records$Training)
    etiology <- unique(records$HF.etiology)
    if (length(arm) != 1L || length(etiology) != 1L) {
      stop("Training and HF.etiology must be constant within patient.",
           call. = FALSE)
    }
    death_records <- records$time[records$status == raw_status["death"]]
    death_time <- if (length(death_records)) min(death_records) else Inf
    censor_records <- records$time[records$status == raw_status["censor"]]
    raw_censor <- if (length(censor_records)) min(censor_records) else tau
    observed_death <- if (death_time <= min(raw_censor, tau)) death_time else Inf
    if (is.finite(observed_death)) {
      visit_limit <- observed_death
      subject_censor <- tau
      censor_event <- FALSE
    } else {
      raw_limit <- min(raw_censor, tau)
      visit_limit <- floor(raw_limit / visit_interval) * visit_interval
      censor_event <- raw_censor < tau
      subject_censor <- if (censor_event) visit_limit else tau
    }
    visits <- seq(0, tau, by = visit_interval)
    visits <- visits[visits <= visit_limit & visits < observed_death]
    if (!length(visits)) {
      visits <- 0
    }
    hospitalizations <- records$time[
      records$status == raw_status["hospitalization"]]
    states <- vapply(visits, function(visit) {
      if (visit == 0) {
        return(unname(ordered_state["no_hospitalization"]))
      }
      if (any(hospitalizations > visit - visit_interval &
              hospitalizations <= visit)) {
        unname(ordered_state["hospitalization"])
      } else {
        unname(ordered_state["no_hospitalization"])
      }
    }, numeric(1))
    subjects[[i]] <- rmtif_subject(
      arm = arm, times = visits, states = states, tau = tau,
      covariates = etiology, censor_time = subject_censor,
      censor_event = censor_event, death_time = observed_death,
      death_state = if (is.finite(observed_death))
        unname(ordered_state["death"]) else NULL,
      id = ids[i]
    )
    subject_summary[[i]] <- data.frame(
      id = ids[i], arm = arm, death = is.finite(observed_death),
      raw_censor_time = raw_censor,
      recurrent_hospitalizations = sum(
        records$status == raw_status["hospitalization"])
    )
  }
  subject_summary <- do.call(rbind, subject_summary)
  arm_summary <- do.call(rbind, lapply(split(subject_summary,
                                             subject_summary$arm), function(x) {
    censor_times <- x$raw_censor_time[!x$death]
    quantiles <- if (length(censor_times)) {
      stats::quantile(censor_times, quantile_probs, names = FALSE)
    } else {
      rep(NA_real_, 3L)
    }
    data.frame(
      arm = x$arm[1L], patients = nrow(x), deaths = sum(x$death),
      non_death_terminal_or_censor_records = sum(!x$death),
      recurrent_hospitalizations = sum(x$recurrent_hospitalizations),
      censor_q1 = quantiles[1L], censor_median = quantiles[2L],
      censor_q3 = quantiles[3L], row.names = NULL
    )
  }))
  rownames(arm_summary) <- NULL
  list(subjects = subjects, subject_summary = subject_summary,
       arm_summary = arm_summary, tau = tau,
       visit_interval = visit_interval)
}

#' Construct scheduled-visit ordered states from recurrent HF event records
#'
#' This entry point accepts the event-record schema used by the public HF
#' example: patient identifier, event time, raw event code, treatment arm, and
#' baseline etiology. Hospitalization since the previous scheduled visit
#' defines the nonterminal state; observed death is an absorbing exact event.
#'
#' @param data Event records with `id`, `time`, `status`, `Training`, and
#'   `HF.etiology` columns.
#' @param tau Restriction horizon.
#' @param visit_interval Scheduled-visit spacing.
#' @param raw_status Named `censor`, `death`, `hospitalization` code vector.
#' @param ordered_state Named ordered-state code vector.
#' @param quantile_probs Three probabilities for the censor-time arm audit.
#' @return List with constructed `rmtif_subject` paths and subject and arm
#'   audit data frames.
#' @export
rmtif_construct_scheduled_visit_states <- function(
    data, tau, visit_interval, raw_status, ordered_state, quantile_probs) {
  rmtif_hfmock_subjects(data, tau, visit_interval, raw_status,
                        ordered_state, quantile_probs)
}

#' Analyze the public `hfmock` scheduled-visit process
#'
#' If `data` is omitted, the function loads `hfmock` only when the optional
#' `Wcompo` package is installed.  No synthetic replacement is silently used.
#' The analysis reports arm counts and terminal-aware IPCW bootstrap inference.
#'
#' @param data Optional `hfmock` data frame.
#' @param tau Restriction horizon.
#' @param visit_interval Scheduled visit spacing.
#' @param bootstrap Stratified bootstrap replicates.
#' @param seed Reproducibility seed.
#' @param raw_status Named raw event-code vector; see `rmtif_hfmock_subjects`.
#' @param ordered_state Named ordered-state vector; see `rmtif_hfmock_subjects`.
#' @param quantile_probs Censor-time audit probabilities; see
#'   `rmtif_hfmock_subjects`.
#' @param conf_level Normal confidence level.
#' @param bootstrap_success_fraction Required fraction of successful bootstrap
#'   resamples before reporting inference.
#' @param minimum_survival Positivity threshold.
#' @param projection_diagnostic Whether to retain the dense fixed-weight
#'   projection diagnostic.
#' @return List with construction audit and an `rmtif_ipcw` fit.
#' @export
rmtif_hfmock_analysis <- function(data = NULL, tau, visit_interval, bootstrap,
                                  seed, raw_status, ordered_state,
                                  quantile_probs, conf_level,
                                  bootstrap_success_fraction,
                                  minimum_survival = 0,
                                  projection_diagnostic = FALSE) {
  if (is.null(data)) {
    if (!requireNamespace("Wcompo", quietly = TRUE)) {
      stop(
        "The optional package 'Wcompo' is required when data is omitted. ",
        "Install it from CRAN or supply the hfmock data frame explicitly.",
        call. = FALSE
      )
    }
    data_environment <- new.env(parent = emptyenv())
    utils::data("hfmock", package = "Wcompo", envir = data_environment)
    if (!exists("hfmock", envir = data_environment, inherits = FALSE)) {
      stop("Wcompo is installed but its hfmock dataset could not be loaded.",
           call. = FALSE)
    }
    data <- get("hfmock", envir = data_environment, inherits = FALSE)
  }
  construction <- rmtif_construct_scheduled_visit_states(
    data, tau, visit_interval, raw_status, ordered_state, quantile_probs)
  fit <- rmtif_ipcw(
    construction$subjects, tau = tau, bootstrap = bootstrap, seed = seed,
    conf_level = conf_level,
    bootstrap_success_fraction = bootstrap_success_fraction,
    minimum_survival = minimum_survival,
    projection_diagnostic = projection_diagnostic
  )
  list(construction = construction, fit = fit)
}


# -----------------------------------------------------------------------------
# Appendix C oncology comparison (version 1.2)
# -----------------------------------------------------------------------------

#' Validate an ordered-state oncology simulation design
#'
#' All baseline intensities, covariate and treatment log-rate effects, state
#' labels, and the initial and terminal states are supplied by the caller.
#' The named clinical states define the response-confirmation endpoints.
#' No manuscript simulation settings are inserted by this constructor.
#'
#' @param transition Data frame with `from`, `to`, `baseline_rate`,
#'   `covariate_log_rate`, and `treatment_log_rate` columns, one row per allowed
#'   transition. Rates are nonnegative; log-rate effects are finite.
#' @param states Named, increasing, zero-based integer codes with names `CR`,
#'   `PR`, `SD`, `PD`, and `death`, in clinical order from better to worse.
#' @param initial_state Initial state code.
#' @param death_state Absorbing terminal state code, equal to `states["death"]`.
#' @param scenario Nonempty scenario label used only in returned metadata.
#' @return A validated list with `scenario`, `transition`, `states`,
#'   `initial_state`, and `death_state`.
#' @export
rmtif_oncology_design <- function(
    transition, states, initial_state, death_state, scenario) {
  state_names <- c("CR", "PR", "SD", "PD", "death")
  if (!is.numeric(states) || !identical(names(states), state_names) ||
      anyNA(states) || any(!is.finite(states)) ||
      !identical(as.integer(states), seq_along(states) - 1L) ||
      any(states != as.integer(states))) {
    stop("states must be increasing zero-based integer codes named CR, PR, SD, PD, death.",
         call. = FALSE)
  }
  initial_state <- .rmtif_scalar_count(initial_state, "initial_state")
  death_state <- .rmtif_scalar_count(death_state, "death_state")
  if (!initial_state %in% states || death_state != unname(states["death"])) {
    stop("initial_state must belong to states and death_state must equal states['death'].",
         call. = FALSE)
  }
  if (!is.character(scenario) || length(scenario) != 1L ||
      is.na(scenario) || !nzchar(scenario)) {
    stop("scenario must be one nonempty character label.", call. = FALSE)
  }
  required <- c("from", "to", "baseline_rate", "covariate_log_rate",
                "treatment_log_rate")
  if (!is.data.frame(transition) || !nrow(transition) ||
      !all(required %in% names(transition)) ||
      !all(vapply(transition[required], is.numeric, logical(1))) ||
      anyNA(transition[required]) ||
      any(!is.finite(as.matrix(transition[required]))) ||
      any(!transition$from %in% states) || any(!transition$to %in% states) ||
      any(transition$from == transition$to) ||
      any(transition$from == death_state) || any(transition$baseline_rate < 0) ||
      anyDuplicated(transition[c("from", "to")])) {
    stop("transition must contain unique allowed nonterminal transitions, finite log-rate effects, and nonnegative baseline rates.",
         call. = FALSE)
  }
  list(
    scenario = scenario, transition = transition,
    states = states, initial_state = initial_state, death_state = death_state
  )
}

# Backward-compatible descriptive alias used by the archived code.
#' @rdname rmtif_oncology_design
#' @export
rmtif_oncology_parameters <- rmtif_oncology_design

#' Calculate conditional oncology transition intensities and jump probabilities
#'
#' Apply the supplied Cox-form transition model at a specified treatment arm
#' and baseline covariate value. The next-state probabilities divide each
#' transition intensity by the total exit rate of its origin state. They are
#' undefined (`NA`) when all allowed exit rates from that state are zero.
#'
#' @param design Validated output of `rmtif_oncology_design()`.
#' @param arm Treatment-arm code, either control `0` or treatment `1`.
#' @param covariate Finite scalar baseline covariate value.
#' @return The design's transition data frame with `intensity`,
#'   `total_exit_rate`, and `next_transition_probability` columns added.
#' @export
rmtif_oncology_transition_rates <- function(design, arm, covariate) {
  design <- do.call(rmtif_oncology_design, design[c(
    "transition", "states", "initial_state", "death_state", "scenario")])
  arm <- .rmtif_scalar_count(arm, "arm")
  if (!arm %in% 0:1) stop("arm must be either zero or one.", call. = FALSE)
  covariate <- .rmtif_scalar_numeric(covariate, "covariate")
  transition <- design$transition
  transition$intensity <- transition$baseline_rate * exp(
    transition$treatment_log_rate * arm +
      transition$covariate_log_rate * covariate)
  if (any(!is.finite(transition$intensity))) {
    stop("Conditional transition intensities overflowed; check rates and effects.",
         call. = FALSE)
  }
  transition$total_exit_rate <- stats::ave(
    transition$intensity, transition$from, FUN = sum)
  transition$next_transition_probability <- ifelse(
    transition$total_exit_rate > 0,
    transition$intensity / transition$total_exit_rate, NA_real_)
  transition
}

.rmtif_oncology_generator <- function(design, arm, covariate) {
  q <- matrix(0, length(design$states), length(design$states),
              dimnames = list(as.character(design$states),
                              as.character(design$states)))
  tr <- design$transition
  value <- tr$baseline_rate * exp(
    tr$treatment_log_rate * arm + tr$covariate_log_rate * covariate
  )
  q[cbind(tr$from + 1L, tr$to + 1L)] <- value
  diag(q) <- -rowSums(q)
  .rmtif_validate_generator(q)
}

.rmtif_oncology_path_state <- function(path, time) {
  # Default findInterval implements a right-continuous path at every jump,
  # including a time exactly equal to the final transition.
  index <- findInterval(time, path$times)
  path$states[pmax(1L, pmin(length(path$states), index))]
}

.rmtif_oncology_endpoints <- function(assessment_times, assessment_states,
                                      death_time, censor_time, states) {
  keep <- assessment_times >= 0 & assessment_times <= censor_time &
    assessment_times < death_time
  assessment_times <- assessment_times[keep]
  assessment_states <- assessment_states[keep]
  post <- which(assessment_times > 0)
  pd <- post[assessment_states[post] == states["PD"]]
  pd_time <- if (length(pd)) assessment_times[pd[1L]] else Inf
  pfs_event_time <- min(pd_time, death_time)
  os_event <- is.finite(death_time) && death_time <= censor_time
  pfs_event <- is.finite(pfs_event_time) && pfs_event_time <= censor_time

  # Confirmation uses consecutive chronological visits, both strictly before
  # PFS and no later than the observation limit.  Baseline SD never confirms
  # disease control.
  eligible <- post[assessment_times[post] < pfs_event_time &
                     assessment_times[post] <= censor_time]
  confirmed_cr <- confirmed_pr <- confirmed_sd <- FALSE
  if (length(eligible) >= 2L) {
    pair_first <- seq_len(length(eligible) - 1L)
    for (k in pair_first) {
      pair <- assessment_states[eligible[c(k, k + 1L)]]
      confirmed_cr <- confirmed_cr || all(pair == states["CR"])
      confirmed_pr <- confirmed_pr || all(pair %in% states[c("CR", "PR")])
      confirmed_sd <- confirmed_sd || all(pair %in% states[c("CR", "PR", "SD")])
    }
  }
  bor <- if (confirmed_cr) 0L else if (confirmed_pr) 1L else
    if (confirmed_sd) 2L else 3L
  list(
    os_time = min(death_time, censor_time), os_event = as.integer(os_event),
    pfs_time = min(pfs_event_time, censor_time),
    pfs_event = as.integer(pfs_event),
    objective_response = as.integer(bor <= 1L),
    bor_ordinal = bor,
    confirmed_disease_control = as.integer(bor <= 2L),
    observed_pd_time = pd_time
  )
}

#' Simulate one shared Appendix C oncology trial
#'
#' All six methods are derived from this same patient-level realization.
#' Nominal visits are generated far enough to cover the jitter window and then
#' visits outside `(0,tau]`, after loss to follow-up, or at/after death are
#' dropped.  The RMT-IF representation is compressed to actual observed state
#' changes; the complete visit history is retained for OR and BOR confirmation.
#'
#' @param n_per_arm Exact subjects in each randomized arm.
#' @param design Validated output of `rmtif_oncology_design()`.
#' @param tau Administrative restriction horizon.
#' @param assessment_interval Nominal imaging interval.
#' @param assessment_jitter Half-width of independent uniform visit jitter.
#' @param censoring_rate Exponential LTFU rate; zero specifies complete follow-up.
#' @param covariate_mean Mean of the Gaussian baseline covariate.
#' @param covariate_sd Positive standard deviation of the Gaussian covariate.
#' @param seed Reproducibility seed.
#' @return Shared endpoint data, RMT-IF subjects, visits, and latent paths.
#' @export
rmtif_simulate_oncology_trial <- function(
    n_per_arm, design, tau, assessment_interval, assessment_jitter,
    censoring_rate, covariate_mean, covariate_sd, seed = NULL) {
  design <- do.call(rmtif_oncology_design, design[c(
    "transition", "states", "initial_state", "death_state", "scenario")])
  n_per_arm <- .rmtif_scalar_count(n_per_arm, "n_per_arm", lower = 2)
  tau <- .rmtif_scalar_numeric(tau, "tau", lower = 0)
  assessment_interval <- .rmtif_scalar_numeric(
    assessment_interval, "assessment_interval", lower = 0)
  assessment_jitter <- .rmtif_scalar_numeric(
    assessment_jitter, "assessment_jitter", lower = 0)
  censoring_rate <- .rmtif_scalar_numeric(
    censoring_rate, "censoring_rate", lower = 0)
  covariate_mean <- .rmtif_scalar_numeric(covariate_mean, "covariate_mean")
  covariate_sd <- .rmtif_scalar_numeric(covariate_sd, "covariate_sd", lower = 0)
  if (covariate_sd <= 0) stop("covariate_sd must be positive.", call. = FALSE)
  if (assessment_interval <= 0 ||
      assessment_jitter >= assessment_interval / 2) {
    stop("assessment_interval must be positive and assessment_jitter must be less than half the interval.",
         call. = FALSE)
  }
  design$tau <- tau
  design$assessment_interval <- assessment_interval
  design$assessment_jitter <- assessment_jitter
  design$censoring_rate <- censoring_rate
  .rmtif_with_seed(seed, {
    arm <- rep(0:1, each = n_per_arm)
    n <- length(arm)
    covariate <- stats::rnorm(n, mean = covariate_mean, sd = covariate_sd)
    loss_time <- if (censoring_rate > 0) {
      stats::rexp(n, censoring_rate)
    } else rep(Inf, n)
    latent <- visits <- subjects <- endpoints <- vector("list", n)
    maximum_nominal <- floor((tau + assessment_jitter) /
                               assessment_interval)
    design$nominal_assessments <- maximum_nominal
    for (i in seq_len(n)) {
      generator <- .rmtif_oncology_generator(design, arm[i], covariate[i])
      latent[[i]] <- rmtif_simulate_ctmc_path(
        generator, design$initial_state, tau)
      death_time <- .rmtif_terminal_time(latent[[i]], design$death_state)
      censor_time <- min(loss_time[i], tau)
      nominal <- assessment_interval * seq_len(maximum_nominal)
      assessment <- nominal + stats::runif(
        maximum_nominal, -assessment_jitter, assessment_jitter)
      # Censoring is evaluated before a tied visit; exact death takes priority.
      assessment <- sort(assessment[
        assessment > 0 & assessment <= tau & assessment < loss_time[i] &
          assessment < death_time
      ])
      assessment_times <- c(0, assessment)
      assessment_states <- .rmtif_oncology_path_state(
        latent[[i]], assessment_times)
      visits[[i]] <- data.frame(
        id = i, time = assessment_times, state = assessment_states)
      endpoints[[i]] <- .rmtif_oncology_endpoints(
        assessment_times, assessment_states, death_time, censor_time,
        states = design$states)

      observed_death <- if (death_time <= censor_time) death_time else Inf
      # Removing repeated observed states is algebraically exact for the
      # piecewise-constant RMT-IF path and greatly reduces its merged grid.
      changed <- c(TRUE, diff(assessment_states) != 0)
      path_times <- assessment_times[changed]
      path_states <- assessment_states[changed]
      subjects[[i]] <- rmtif_subject(
        arm = arm[i], times = path_times, states = path_states, tau = tau,
        covariates = covariate[i], censor_time = censor_time,
        censor_event = loss_time[i] < min(death_time, tau),
        death_time = observed_death,
        death_state = if (is.finite(observed_death)) design$death_state else NULL,
        id = i)
    }
    endpoint <- data.frame(
      id = seq_len(n), arm = arm, covariate = covariate,
      os_time = vapply(endpoints, `[[`, numeric(1), "os_time"),
      os_event = vapply(endpoints, `[[`, integer(1), "os_event"),
      pfs_time = vapply(endpoints, `[[`, numeric(1), "pfs_time"),
      pfs_event = vapply(endpoints, `[[`, integer(1), "pfs_event"),
      objective_response = vapply(
        endpoints, `[[`, integer(1), "objective_response"),
      bor_ordinal = vapply(endpoints, `[[`, integer(1), "bor_ordinal"),
      confirmed_disease_control = vapply(
        endpoints, `[[`, integer(1), "confirmed_disease_control"),
      ltfu = as.integer(loss_time < pmin(
        vapply(latent, .rmtif_terminal_time, numeric(1),
               death_state = design$death_state), tau)),
      stringsAsFactors = FALSE)
    list(subjects = subjects, data = endpoint, visits = visits,
         latent = latent, design = design, scenario = design$scenario, tau = tau,
         seed = seed)
  })
}

.rmtif_oncology_result <- function(method, estimand, estimate = NA_real_,
                                    standard_error = NA_real_,
                                    statistic = NA_real_, p_value = NA_real_,
                                    converged = FALSE, failure_reason = "",
                                    warning = "", inference = "model Wald",
                                    bootstrap_attempted = 0L,
                                    bootstrap_success = 0L) {
  data.frame(
    method = method, estimand = estimand, estimate = as.numeric(estimate),
    standard_error = as.numeric(standard_error),
    statistic = as.numeric(statistic), p_value = as.numeric(p_value),
    converged = as.logical(converged),
    failure_reason = as.character(failure_reason),
    warning = as.character(warning), inference = as.character(inference),
    bootstrap_attempted = as.integer(bootstrap_attempted),
    bootstrap_success = as.integer(bootstrap_success),
    stringsAsFactors = FALSE)
}

.rmtif_oncology_warnings <- function(code) {
  messages <- character()
  value <- withCallingHandlers(
    tryCatch(code, error = identity),
    warning = function(w) {
      messages <<- c(messages, conditionMessage(w))
      invokeRestart("muffleWarning")
    })
  list(value = value, warning = paste(unique(messages), collapse = " | "))
}

.rmtif_oncology_cox <- function(data, time, event, method) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    return(.rmtif_oncology_result(
      method, "log hazard ratio for treatment", failure_reason =
        "Package 'survival' is not installed."))
  }
  formula <- stats::as.formula(paste0(
    "survival::Surv(", time, ",", event, ") ~ arm + covariate"))
  captured <- .rmtif_oncology_warnings(survival::coxph(
    formula, data = data, ties = "efron", singular.ok = FALSE,
    model = FALSE, x = FALSE, y = FALSE))
  if (inherits(captured$value, "error")) {
    return(.rmtif_oncology_result(
      method, "log hazard ratio for treatment",
      failure_reason = conditionMessage(captured$value),
      warning = captured$warning))
  }
  fit <- captured$value
  estimate <- unname(stats::coef(fit)["arm"])
  standard_error <- tryCatch(
    sqrt(stats::vcov(fit)["arm", "arm"]), error = function(e) NA_real_)
  failed_warning <- grepl("converg|infinite|singular|ran out", captured$warning,
                         ignore.case = TRUE)
  good <- is.null(fit$fail) && !failed_warning && is.finite(estimate) &&
    is.finite(standard_error) && standard_error > 0
  statistic <- if (good) estimate / standard_error else NA_real_
  .rmtif_oncology_result(
    method, "log hazard ratio for treatment", estimate, standard_error,
    statistic, if (good) 2 * stats::pnorm(-abs(statistic)) else NA_real_,
    converged = good,
    failure_reason = if (good) "" else
      "Cox fit failed convergence checks or returned a non-finite coefficient/SE.",
    warning = captured$warning,
    inference = "Cox model-based Wald; Efron ties")
}

.rmtif_oncology_logistic <- function(data, response, method, estimand) {
  formula <- stats::as.formula(paste(response, "~ arm + covariate"))
  captured <- .rmtif_oncology_warnings(stats::glm(
    formula, data = data, family = stats::binomial()))
  if (inherits(captured$value, "error")) {
    return(.rmtif_oncology_result(
      method, estimand, failure_reason = conditionMessage(captured$value),
      warning = captured$warning))
  }
  fit <- captured$value
  coefficient <- summary(fit)$coefficients
  estimate <- if ("arm" %in% rownames(coefficient)) coefficient["arm", 1L] else NA_real_
  standard_error <- if ("arm" %in% rownames(coefficient)) coefficient["arm", 2L] else NA_real_
  full_rank <- fit$rank == length(stats::coef(fit))
  separated <- length(unique(data[[response]])) < 2L ||
    any(stats::fitted(fit) <= sqrt(.Machine$double.eps) |
          stats::fitted(fit) >= 1 - sqrt(.Machine$double.eps)) ||
    grepl("fitted probabilities numerically 0 or 1|did not converge",
          captured$warning, ignore.case = TRUE)
  good <- isTRUE(fit$converged) && !separated && full_rank && is.finite(estimate) &&
    is.finite(standard_error) && standard_error > 0
  statistic <- if (good) estimate / standard_error else NA_real_
  .rmtif_oncology_result(
    method, estimand, estimate, standard_error, statistic,
    if (good) 2 * stats::pnorm(-abs(statistic)) else NA_real_, good,
    failure_reason = if (good) "" else
      "Logistic fit failed convergence/rank checks or showed degeneracy or near-boundary fitted probabilities.",
    warning = captured$warning,
    inference = "binomial-logit model-based Wald")
}

.rmtif_oncology_bor <- function(data) {
  estimand <- paste(
    "study-specific proportional-odds treatment coefficient;",
    "positive MASS::polr coefficient shifts toward worse BOR")
  if (!requireNamespace("MASS", quietly = TRUE)) {
    return(.rmtif_oncology_result(
      "BOR", estimand, failure_reason = "Package 'MASS' is not installed."))
  }
  observed <- sort(unique(data$bor_ordinal))
  if (length(observed) < 3L) {
    return(.rmtif_oncology_result(
      "BOR", estimand,
      failure_reason = "Fewer than three ordinal BOR categories were observed."))
  }
  model_data <- data
  model_data$bor <- ordered(model_data$bor_ordinal, levels = 0:3)
  captured <- .rmtif_oncology_warnings(MASS::polr(
    bor ~ arm + covariate, data = model_data, method = "logistic",
    Hess = TRUE, model = FALSE))
  if (inherits(captured$value, "error")) {
    return(.rmtif_oncology_result(
      "BOR", estimand, failure_reason = conditionMessage(captured$value),
      warning = captured$warning))
  }
  fit <- captured$value
  coefficient <- tryCatch(summary(fit)$coefficients,
                          error = function(e) NULL)
  estimate <- if (!is.null(coefficient) && "arm" %in% rownames(coefficient)) {
    coefficient["arm", "Value"]
  } else NA_real_
  standard_error <- if (!is.null(coefficient) &&
                        "arm" %in% rownames(coefficient)) {
    coefficient["arm", "Std. Error"]
  } else NA_real_
  good <- identical(fit$convergence, 0L) && is.finite(estimate) &&
    is.finite(standard_error) && standard_error > 0
  statistic <- if (good) estimate / standard_error else NA_real_
  .rmtif_oncology_result(
    "BOR", estimand, estimate, standard_error, statistic,
    if (good) 2 * stats::pnorm(-abs(statistic)) else NA_real_, good,
    failure_reason = if (good) "" else
      "Proportional-odds fit failed or had a non-finite treatment SE.",
    warning = captured$warning,
    inference = "MASS::polr proportional-odds model-based Wald")
}

.rmtif_oncology_rmt <- function(simulation,
                                 inference = c("complete_projection",
                                               "ipcw_bootstrap"),
                                 bootstrap, bootstrap_seed, conf_level,
                                 bootstrap_success_fraction,
                                 minimum_survival) {
  inference <- match.arg(inference)
  censoring <- if (inference == "complete_projection") "none" else "arm_km"
  bootstrap_requested <- if (inference == "complete_projection") 0L else bootstrap
  captured <- .rmtif_oncology_warnings(rmtif_regression(
    simulation$subjects, tau = simulation$tau, censoring = censoring,
    bootstrap = bootstrap_requested, conf_level = conf_level,
    bootstrap_success_fraction = bootstrap_success_fraction,
    seed = bootstrap_seed, minimum_survival = minimum_survival))
  if (inherits(captured$value, "error")) {
    return(.rmtif_oncology_result(
      "RMT-IF", "pairwise treatment coefficient at equal covariates",
      failure_reason = conditionMessage(captured$value),
      warning = captured$warning, bootstrap_attempted = bootstrap_requested))
  }
  fit <- captured$value
  estimate <- unname(fit$coefficients["treatment"])
  treatment_index <- match("treatment", names(fit$coefficients))
  standard_error <- unname(fit$se[treatment_index])
  statistic <- unname(fit$statistic["treatment"])
  p_value <- unname(fit$p_value["treatment"])
  good <- isTRUE(fit$root_converged) && all(is.finite(c(
    estimate, standard_error, statistic, p_value))) && standard_error > 0
  .rmtif_oncology_result(
    "RMT-IF", "pairwise treatment coefficient at equal covariates",
    estimate, standard_error, statistic, p_value, good,
    failure_reason = if (good) "" else
      "RMT-IF did not return a verified root and finite treatment inference.",
    warning = captured$warning, inference = fit$inference,
    bootstrap_attempted = bootstrap_requested,
    bootstrap_success = fit$bootstrap_success)
}

#' Compare six covariate-adjusted tests in one oncology trial
#'
#' OS/PFS use Cox Wald tests, OR and CDC use logistic Wald tests, study-specific
#' ordinal BOR uses a proportional-odds Wald test, and RMT-IF tests its pairwise
#' treatment coefficient.  Every method uses the same generated subjects.
#'
#' @param simulation Output of `rmtif_simulate_oncology_trial()`.
#' @param rmt_inference Complete-path projection or censored IPCW bootstrap.
#' @param regression_bootstrap IPCW bootstrap count (ignored for complete data).
#' @param bootstrap_seed Separate seed for the stratified bootstrap.
#' @param conf_level Confidence level.
#' @param bootstrap_success_fraction Required successful bootstrap fraction.
#' @param minimum_survival IPCW positivity threshold.
#' @param test_level Two-sided rejection threshold.
#' @return Six-row long-format method result.
#' @export
rmtif_oncology_compare_trial <- function(
    simulation, rmt_inference, regression_bootstrap,
    conf_level, bootstrap_success_fraction, minimum_survival, test_level,
    bootstrap_seed = NULL) {
  rmt_inference <- match.arg(rmt_inference,
                             c("complete_projection", "ipcw_bootstrap"))
  test_level <- .rmtif_conf_level(test_level)
  data <- simulation$data
  if (rmt_inference == "complete_projection" && any(data$ltfu > 0L)) {
    stop("Complete-path projection requires no loss to follow-up; use ipcw_bootstrap for censored trials.",
         call. = FALSE)
  }
  results <- rbind(
    .rmtif_oncology_cox(data, "os_time", "os_event", "OS"),
    .rmtif_oncology_cox(data, "pfs_time", "pfs_event", "PFS"),
    .rmtif_oncology_logistic(
      data, "objective_response", "OR",
      "log odds ratio for confirmed objective response"),
    .rmtif_oncology_bor(data),
    .rmtif_oncology_logistic(
      data, "confirmed_disease_control", "CDC",
      "log odds ratio for confirmed disease control"),
    .rmtif_oncology_rmt(
      simulation, rmt_inference, regression_bootstrap, bootstrap_seed,
      conf_level, bootstrap_success_fraction, minimum_survival))
  results$reject <- results$converged & is.finite(results$p_value) &
    results$p_value < test_level
  results
}

#' Monte Carlo comparison of oncology endpoints on shared trials
#'
#' Simulate each scenario using a separate seed for every trial and bootstrap.
#' All six methods analyse the same subjects within a replicate. Failed fits
#' remain in the replicate file; rejection proportions and their Monte Carlo
#' standard errors use the number of valid fits as the denominator. Optional
#' scenario checkpoints permit an interrupted run to resume without changing
#' the random streams. Parallel workers are allocated by scenario, not method.
#'
#' @param scenarios Distinct scenario labels to simulate.
#' @param designs Named list of validated `rmtif_oncology_design()` results,
#'   with one entry per requested scenario. List names and each design's
#'   `scenario` label must agree.
#' @param n_per_arm Number of subjects in each arm.
#' @param replications Number of attempted trials per scenario.
#' @param tau Administrative horizon, in months.
#' @param assessment_interval Nominal interval between assessments.
#' @param assessment_jitter Half-width of uniform assessment-time jitter.
#' @param censoring_rate Exponential loss-to-follow-up rate; zero gives complete
#'   follow-up through the administrative horizon.
#' @param covariate_mean Mean of the Gaussian baseline covariate.
#' @param covariate_sd Positive standard deviation of the Gaussian covariate.
#' @param rmt_inference Complete-path projection or censored IPCW bootstrap.
#' @param regression_bootstrap Bootstrap resamples for censored RMT-IF inference.
#' @param seed Master seed; replicate seeds are retained in the output.
#' @param conf_level Confidence level for fitted models.
#' @param bootstrap_success_fraction Required successful bootstrap fraction.
#' @param minimum_survival IPCW positivity threshold.
#' @param test_level Two-sided rejection threshold.
#' @param workers Number of scenario-level PSOCK workers.
#' @param checkpoint_dir Optional directory for resumable scenario checkpoints.
#' @return List containing summary, replicate, transition-design, and runtime
#'   data frames. The example driver saves these outside the package.
#' @export
rmtif_oncology_study <- function(
    scenarios, designs, n_per_arm, replications, tau, assessment_interval,
    assessment_jitter, censoring_rate, covariate_mean, covariate_sd,
    rmt_inference, regression_bootstrap, seed, conf_level,
    bootstrap_success_fraction, minimum_survival, test_level,
    workers = 1L, checkpoint_dir = NULL) {
  rmt_inference <- match.arg(rmt_inference,
                             c("complete_projection", "ipcw_bootstrap"))
  if (!is.character(scenarios) || !length(scenarios) || anyNA(scenarios) ||
      any(!nzchar(scenarios)) || anyDuplicated(scenarios) ||
      !is.list(designs) || anyDuplicated(names(designs)) ||
      !all(scenarios %in% names(designs))) {
    stop("scenarios must be distinct labels with matching named entries in designs.",
         call. = FALSE)
  }
  designs <- lapply(designs[scenarios], function(design) {
    do.call(rmtif_oncology_design, design[c(
      "transition", "states", "initial_state", "death_state", "scenario")])
  })
  if (!identical(vapply(designs, `[[`, character(1), "scenario"),
                 stats::setNames(scenarios, scenarios))) {
    stop("Each design's scenario label must match its list name.", call. = FALSE)
  }
  n_per_arm <- .rmtif_scalar_count(n_per_arm, "n_per_arm", lower = 2)
  tau <- .rmtif_scalar_numeric(tau, "tau", lower = 0)
  assessment_interval <- .rmtif_scalar_numeric(
    assessment_interval, "assessment_interval", lower = 0)
  assessment_jitter <- .rmtif_scalar_numeric(
    assessment_jitter, "assessment_jitter", lower = 0)
  if (assessment_interval <= 0 ||
      assessment_jitter >= assessment_interval / 2) {
    stop("assessment_interval must be positive and assessment_jitter must be less than half the interval.",
         call. = FALSE)
  }
  censoring_rate <- .rmtif_scalar_numeric(censoring_rate, "censoring_rate", lower = 0)
  covariate_mean <- .rmtif_scalar_numeric(covariate_mean, "covariate_mean")
  covariate_sd <- .rmtif_scalar_numeric(covariate_sd, "covariate_sd", lower = 0)
  if (covariate_sd <= 0) stop("covariate_sd must be positive.", call. = FALSE)
  conf_level <- .rmtif_conf_level(conf_level)
  bootstrap_success_fraction <- .rmtif_bootstrap_success_fraction(
    bootstrap_success_fraction)
  minimum_survival <- .rmtif_minimum_survival(minimum_survival)
  regression_bootstrap <- .rmtif_scalar_count(
    regression_bootstrap, "regression_bootstrap")
  if (rmt_inference == "ipcw_bootstrap" && regression_bootstrap < 2L) {
    stop("ipcw_bootstrap inference requires at least two bootstrap resamples.",
         call. = FALSE)
  }
  if (length(replications) != 1L || !is.finite(replications) ||
      replications < 1 || replications != as.integer(replications)) {
    stop("replications must be a positive integer.")
  }
  if (length(workers) != 1L || !is.finite(workers) || workers < 1 ||
      workers != as.integer(workers)) stop("workers must be a positive integer.")
  if (!is.finite(test_level) || test_level <= 0 || test_level >= 1) {
    stop("test_level must lie strictly between zero and one.")
  }
  if (rmt_inference == "complete_projection" && censoring_rate != 0) {
    stop("Complete-path inference requires censoring_rate = 0.")
  }
  replications <- as.integer(replications)
  workers <- min(as.integer(workers), length(scenarios))
  methods <- c("OS", "PFS", "OR", "BOR", "CDC", "RMT-IF")
  seeds <- .rmtif_with_seed(seed, matrix(sample.int(
    .Machine$integer.max, 2L * length(scenarios) * replications), ncol = 2L))
  signature <- list(scenarios = scenarios, designs = designs, n_per_arm = n_per_arm,
    replications = replications, tau = tau,
    assessment_interval = assessment_interval,
    assessment_jitter = assessment_jitter, censoring_rate = censoring_rate,
    covariate_mean = covariate_mean, covariate_sd = covariate_sd,
    rmt_inference = rmt_inference, regression_bootstrap = regression_bootstrap,
    seed = seed, conf_level = conf_level,
    bootstrap_success_fraction = bootstrap_success_fraction,
    minimum_survival = minimum_survival, test_level = test_level,
    seeds = seeds, implementation = "oncology-v1.2-explicit-design-admin-inclusive")
  if (!is.null(checkpoint_dir)) {
    dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
    checkpoint_dir <- normalizePath(checkpoint_dir, mustWork = TRUE)
  }
  started <- Sys.time()
  run_scenario <- function(scenario_index) {
    scenario <- scenarios[scenario_index]
    rows <- vector("list", replications)
    checkpoint <- if (is.null(checkpoint_dir)) NULL else file.path(
      checkpoint_dir, paste0("oncology_", scenario, "_v1.2.rds"))
    completed <- 0L
    if (!is.null(checkpoint) && file.exists(checkpoint)) {
      previous <- readRDS(checkpoint)
      if (!identical(previous$signature, signature)) {
        stop("Checkpoint settings differ: ", checkpoint)
      }
      completed <- previous$completed
      rows <- previous$rows
    }
    scenario_start <- Sys.time()
    if (completed < replications) for (replicate in seq.int(
        completed + 1L, replications)) {
      seed_row <- (scenario_index - 1L) * replications + replicate
      outer_seed <- seeds[seed_row, 1L]
      bootstrap_seed <- seeds[seed_row, 2L]
      trial_start <- Sys.time()
      trial <- tryCatch(rmtif_simulate_oncology_trial(
        n_per_arm = n_per_arm, design = designs[[scenario]], tau = tau,
        assessment_interval = assessment_interval,
        assessment_jitter = assessment_jitter,
        censoring_rate = censoring_rate, covariate_mean = covariate_mean,
        covariate_sd = covariate_sd, seed = outer_seed), error = identity)
      result <- tryCatch({
        if (inherits(trial, "error")) stop(conditionMessage(trial))
        rmtif_oncology_compare_trial(trial, rmt_inference = rmt_inference,
          regression_bootstrap = regression_bootstrap,
          bootstrap_seed = bootstrap_seed, conf_level = conf_level,
          bootstrap_success_fraction = bootstrap_success_fraction,
          minimum_survival = minimum_survival, test_level = test_level)
      }, error = function(e) {
        failed <- do.call(rbind, lapply(methods, function(method) {
          .rmtif_oncology_result(method, "unavailable",
            failure_reason = conditionMessage(e))
        }))
        failed$reject <- FALSE
        failed
      })
      result$scenario <- scenario
      result$replicate <- replicate
      result$outer_seed <- outer_seed
      result$bootstrap_seed <- bootstrap_seed
      result$n_control <- n_per_arm
      result$n_treatment <- n_per_arm
      result$tau <- tau
      for (arm in 0:1) {
        result[[paste0("ltfu_arm", arm)]] <- if (inherits(trial, "error")) {
          NA_integer_
        } else sum(trial$data$ltfu[trial$data$arm == arm])
        for (category in 0:3) {
          result[[paste0("bor", category, "_arm", arm)]] <-
            if (inherits(trial, "error")) NA_integer_ else sum(
              trial$data$bor_ordinal[trial$data$arm == arm] == category)
        }
      }
      result$trial_elapsed_seconds <- as.numeric(difftime(
        Sys.time(), trial_start, units = "secs"))
      rows[[replicate]] <- result
      if (replicate %% 25L == 0L || replicate == replications) {
        if (!is.null(checkpoint)) saveRDS(list(signature = signature,
          completed = replicate, rows = rows), checkpoint)
      }
    }
    list(replicates = do.call(rbind, rows), runtime = data.frame(
      scenario = scenario, resumed_after = completed,
      elapsed_seconds = as.numeric(difftime(Sys.time(), scenario_start,
                                           units = "secs"))))
  }
  if (workers == 1L) {
    runs <- lapply(seq_along(scenarios), run_scenario)
  } else {
    cluster <- parallel::makePSOCKcluster(workers)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    function_environment <- environment(rmtif_oncology_study)
    exports <- ls(function_environment, all.names = TRUE)
    exports <- exports[grepl("^\\.?rmtif_", exports)]
    parallel::clusterExport(cluster, exports, envir = function_environment)
    runs <- parallel::parLapply(cluster, seq_along(scenarios), run_scenario)
  }
  replicates <- do.call(rbind, lapply(runs, `[[`, "replicates"))
  rownames(replicates) <- NULL
  summary <- do.call(rbind, lapply(scenarios, function(scenario) {
    do.call(rbind, lapply(methods, function(method) {
      subset <- replicates[replicates$scenario == scenario &
                             replicates$method == method, , drop = FALSE]
      valid <- subset$converged & is.finite(subset$p_value)
      successful <- sum(valid)
      rejections <- sum(subset$reject[valid])
      proportion <- if (successful) rejections / successful else NA_real_
      reasons <- unique(subset$failure_reason[!valid])
      data.frame(scenario = scenario, method = method,
        attempted = nrow(subset), successful = successful,
        failed = sum(!valid), rejections = rejections,
        rejection_proportion = proportion,
        monte_carlo_se = sqrt(proportion * (1 - proportion) / successful),
        rejection_per_attempt = rejections / nrow(subset),
        mean_estimate = if (successful) mean(subset$estimate[valid]) else NA_real_,
        mean_standard_error = if (successful) mean(subset$standard_error[valid])
          else NA_real_, warning_trials = sum(nzchar(subset$warning)),
        failure_reasons = paste(reasons[nzchar(reasons)], collapse = " | "))
    }))
  }))
  design <- do.call(rbind, lapply(scenarios, function(scenario) {
    transition <- designs[[scenario]]$transition
    transition$scenario <- scenario
    transition$treatment_multiplier <- exp(transition$treatment_log_rate)
    transition$n_per_arm <- n_per_arm
    transition$replications <- replications
    transition$tau <- tau
    transition$assessment_interval <- assessment_interval
    transition$assessment_jitter <- assessment_jitter
    transition$censoring_rate <- censoring_rate
    transition$covariate_mean <- covariate_mean
    transition$covariate_sd <- covariate_sd
    transition$initial_state <- designs[[scenario]]$initial_state
    transition$death_state <- designs[[scenario]]$death_state
    transition$rmt_inference <- rmt_inference
    transition$regression_bootstrap <- if (rmt_inference ==
      "complete_projection") 0L else regression_bootstrap
    transition$test_level <- test_level
    transition$seed <- seed
    transition$R_version <- R.version.string
    transition$RNG_kind <- paste(RNGkind(), collapse = "; ")
    transition
  }))
  runtime <- do.call(rbind, lapply(runs, `[[`, "runtime"))
  runtime <- rbind(runtime, data.frame(scenario = "Total wall time",
    resumed_after = NA_integer_, elapsed_seconds = as.numeric(difftime(
      Sys.time(), started, units = "secs"))))
  list(summary = summary, replicates = replicates,
       design = design, runtime = runtime)
}


