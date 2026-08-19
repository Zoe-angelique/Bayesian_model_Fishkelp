# Functions for Zoë's habitat and trophic evidence model

integrate_question_priors <- function(question_answers) {
  alpha <- sum(question_answers$weight * question_answers$alpha)
  beta <- sum(question_answers$weight * question_answers$beta)

  data.frame(
    alpha = alpha,
    beta = beta,
    mean = alpha / (alpha + beta),
    lower_95 = qbeta(0.025, alpha, beta),
    upper_95 = qbeta(0.975, alpha, beta)
  )
}

draw_integrated_priors <- function(
  question_answers,
  weight_concentration,
  n_draws,
  seed
) {
  set.seed(seed)

  n_questions <- nrow(question_answers)
  dirichlet_shape <- weight_concentration * question_answers$weight
  raw_weights <- matrix(
    rgamma(
      n_draws * n_questions,
      shape = rep(dirichlet_shape, times = n_draws)
    ),
    nrow = n_draws,
    byrow = TRUE
  )
  sampled_weights <- raw_weights / rowSums(raw_weights)

  data.frame(
    alpha = as.vector(sampled_weights %*% question_answers$alpha),
    beta = as.vector(sampled_weights %*% question_answers$beta)
  )
}

summarise_integrated_prior_draws <- function(prior_draws, route) {
  prior_mean <- prior_draws$alpha / (prior_draws$alpha + prior_draws$beta)

  data.frame(
    route = route,
    mean_rank = mean(prior_mean),
    lower_95 = unname(quantile(prior_mean, 0.025)),
    upper_95 = unname(quantile(prior_mean, 0.975))
  )
}

fit_evidence_route <- function(
  records,
  prior_draws,
  specificity_sd,
  source_variation_scale = 0.20
) {
  theta <- seq(0.0005, 0.9995, length.out = 2500)
  theta_width <- theta[2] - theta[1]
  source_variation <- seq(0, 0.60, length.out = 121)
  source_variation_width <- source_variation[2] - source_variation[1]

  # These boundaries divide the continuous HDR or TDR into six evidence levels.
  boundaries <- c(-Inf, 0.10, 0.30, 0.50, 0.70, 0.90, Inf)
  lower_bound <- boundaries[records$level + 1]
  upper_bound <- boundaries[records$level + 2]
  specificity_error <- unname(specificity_sd[records$specificity])

  # Average across plausible question weights to obtain the marginal prior.
  prior_density <- vapply(
    theta,
    function(value) {
      mean(dbeta(value, prior_draws$alpha, prior_draws$beta))
    },
    numeric(1)
  )

  # The source-variation parameter measures agreement among literature records.
  # Agreement supports a small value; disagreement supports a larger value.
  log_joint_posterior <- vapply(
    source_variation,
    function(variation) {
      observation_sd <- sqrt(specificity_error^2 + variation^2)

      log_likelihood <- vapply(
        theta,
        function(value) {
          record_probability <- pnorm(
            (upper_bound - value) / observation_sd
          ) - pnorm(
            (lower_bound - value) / observation_sd
          )

          sum(log(pmax(record_probability, .Machine$double.xmin)))
        },
        numeric(1)
      )

      log(prior_density) +
        log(2) +
        dnorm(
          variation,
          mean = 0,
          sd = source_variation_scale,
          log = TRUE
        ) +
        log_likelihood
    },
    numeric(length(theta))
  )

  joint_density <- exp(log_joint_posterior - max(log_joint_posterior))
  joint_density <- joint_density /
    sum(joint_density * theta_width * source_variation_width)

  posterior_density <- rowSums(joint_density) * source_variation_width
  source_variation_density <- colSums(joint_density) * theta_width

  list(
    theta = theta,
    prior_density = prior_density,
    posterior_density = posterior_density,
    grid_width = theta_width,
    source_variation = source_variation,
    source_variation_density = source_variation_density,
    source_variation_width = source_variation_width
  )
}

summarise_evidence_route <- function(fit, route) {
  posterior_mass <- fit$posterior_density * fit$grid_width
  posterior_cdf <- cumsum(posterior_mass)
  credible_limits <- approx(
    x = posterior_cdf,
    y = fit$theta,
    xout = c(0.025, 0.50, 0.975),
    ties = "ordered",
    rule = 2
  )$y

  source_variation_mass <- fit$source_variation_density *
    fit$source_variation_width
  source_variation_cdf <- cumsum(source_variation_mass)
  source_variation_limits <- approx(
    x = source_variation_cdf,
    y = fit$source_variation,
    xout = c(0.025, 0.50, 0.975),
    ties = "ordered",
    rule = 2
  )$y

  data.frame(
    route = route,
    posterior_mean = sum(fit$theta * posterior_mass),
    posterior_median = credible_limits[2],
    lower_95 = credible_limits[1],
    upper_95 = credible_limits[3],
    expected_level = 5 * sum(fit$theta * posterior_mass),
    expected_level_lower_95 = 5 * credible_limits[1],
    expected_level_upper_95 = 5 * credible_limits[3],
    source_variation_median = source_variation_limits[2],
    source_variation_lower_95 = source_variation_limits[1],
    source_variation_upper_95 = source_variation_limits[3]
  )
}

calculate_level_probabilities <- function(fit, route) {
  level_number <- cut(
    fit$theta,
    breaks = c(0, 0.10, 0.30, 0.50, 0.70, 0.90, 1),
    labels = 0:5,
    include.lowest = TRUE,
    right = FALSE
  )

  level_names <- c(
    "None",
    "Co-occurrence",
    "Use",
    "Selection",
    "Contribution",
    "Dependence"
  )

  level_mass <- aggregate(
    fit$posterior_density * fit$grid_width,
    by = list(level = level_number),
    FUN = sum
  )

  data.frame(
    route = route,
    level = 0:5,
    evidence_level = level_names,
    posterior_probability = level_mass$x
  )
}
