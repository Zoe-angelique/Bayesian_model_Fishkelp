#!/usr/bin/env Rscript

# Refit the FishKelp v7 model and extract posterior evidentiary-reach draws for
# every modelled fish listed in Blamey & Bolton (2018), Table 2.
#
# Run from the repository root:
#
#   Rscript R/run_fishkelp_v7_suite.R
#
# IMPORTANT: the hierarchical model is fitted to the complete FishKelp model
# object. The output is filtered to the Blamey & Bolton fish only after fitting.

repo_root <- normalizePath(".", mustWork = TRUE)
bundle_dir <- file.path(repo_root, "data", "fishkelp_v7")
model_data_path <- file.path(bundle_dir, "model", "tdr_model_data_v7.rds")
fish_path <- file.path(bundle_dir, "blamey_bolton_2018_fish.csv")
authoritative_path <- file.path(
  bundle_dir, "fishkelp_posterior_estimates_v7.csv"
)
sampler_path <- file.path(repo_root, "R", "tdr_sampler_v7.R")
output_dir <- file.path(repo_root, "model_outputs")
figure_data_dir <- file.path(repo_root, "Bayesian_framework", "posterior_outputs")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figure_data_dir, recursive = TRUE, showWarnings = FALSE)

source(sampler_path)
model_data <- readRDS(model_data_path)
paper_fish <- read.csv(fish_path, stringsAsFactors = FALSE, check.names = FALSE)
authoritative <- read.csv(
  authoritative_path,
  stringsAsFactors = FALSE,
  check.names = FALSE,
  na.strings = c("NA", "")
)

target_species <- paper_fish$accepted_name[paper_fish$modelled_v7]
target_indices <- match(target_species, model_data$species_lab)

cat(sprintf(
  paste0(
    "Fitting FishKelp v7 to %d records and %d species; ",
    "posterior output will retain %d Blamey & Bolton fish.\n"
  ),
  model_data$N,
  model_data$I,
  length(target_species)
))
cat("Production settings: 4 chains x 60,000 iterations; burn-in 15,000.\n")

chains <- run_model(
  dat = model_data,
  y = model_data$y,
  dimension = model_data$dimension,
  species = model_data$species,
  sigma = model_data$sigma,
  X = model_data$X,
  nchains = 4,
  niter = 60000,
  nburn = 15000,
  nu0 = 40,
  tau0 = 0.7,
  cut_sd = 5
)

# Check every free fixed effect and free cut-point.
all_rhat <- unlist(lapply(seq_len(model_data$D), function(route_index) {
  c(
    sapply(seq_len(model_data$P), function(predictor_index) {
      draws <- sapply(chains, function(chain) {
        chain$beta[, route_index, predictor_index]
      })
      mcmc_diag(draws)[["Rhat"]]
    }),
    sapply(2:(model_data$C - 1L), function(cutpoint_index) {
      draws <- sapply(chains, function(chain) {
        chain$kappa[, route_index, cutpoint_index]
      })
      mcmc_diag(draws)[["Rhat"]]
    })
  )
}))

max_rhat <- max(all_rhat)
if (max_rhat > 1.05) {
  stop(sprintf("Convergence failed: maximum split-Rhat = %.3f", max_rhat))
}

write.csv(
  data.frame(
    n_free_parameters = length(all_rhat),
    max_split_rhat = max_rhat
  ),
  file.path(output_dir, "tdr_convergence_rerun_v7.csv"),
  row.names = FALSE
)

reach_from_chain <- function(chain, species_index, route_index) {
  location <- as.vector(
    chain$beta[, route_index, ] %*% model_data$X_ref[species_index,]
  ) + chain$u[, species_index, route_index]
  cutpoints <- matrix(
    chain$kappa[, route_index,],
    nrow = length(location)
  )

  rowSums(1 - pnorm(cutpoints - location))
}

n_saved_per_chain <- dim(chains[[1]]$beta)[1]
n_draws <- n_saved_per_chain * length(chains)
reach_draws <- array(
  NA_real_,
  dim = c(n_draws, length(target_species), model_data$D),
  dimnames = list(
    NULL,
    target_species,
    c("Trophic", "Habitat")
  )
)

for (chain_index in seq_along(chains)) {
  draw_rows <- ((chain_index - 1L) * n_saved_per_chain + 1L):
    (chain_index * n_saved_per_chain)
  for (species_position in seq_along(target_indices)) {
    for (route_index in seq_len(model_data$D)) {
      reach_draws[draw_rows, species_position, route_index] <- reach_from_chain(
        chains[[chain_index]],
        target_indices[[species_position]],
        route_index
      )
    }
  }
}

summarise_draws <- function(values) {
  c(
    mean = mean(values),
    median = median(values),
    q2.5 = unname(quantile(values, 0.025)),
    q97.5 = unname(quantile(values, 0.975))
  )
}

summary_rows <- do.call(rbind, lapply(seq_along(target_species), function(i) {
  trophic <- summarise_draws(reach_draws[, i, "Trophic"])
  habitat <- summarise_draws(reach_draws[, i, "Habitat"])
  data.frame(
    species = target_species[[i]],
    T_reach = trophic[["mean"]],
    T_median = trophic[["median"]],
    T_reach_lo = trophic[["q2.5"]],
    T_reach_hi = trophic[["q97.5"]],
    H_reach = habitat[["mean"]],
    H_median = habitat[["median"]],
    H_reach_lo = habitat[["q2.5"]],
    H_reach_hi = habitat[["q97.5"]],
    stringsAsFactors = FALSE
  )
}))

summary_rows <- summary_rows[order(-summary_rows$T_reach),]
write.csv(
  summary_rows,
  file.path(output_dir, "blamey_bolton_fish_posterior_summary_v7.csv"),
  row.names = FALSE
)

authoritative_modelled <- authoritative[authoritative$modelled_v7,]
authoritative_modelled <- authoritative_modelled[
  match(summary_rows$species, authoritative_modelled$accepted_name),
]
comparison_columns <- c(
  "T_reach", "T_reach_lo", "T_reach_hi",
  "H_reach", "H_reach_lo", "H_reach_hi"
)
largest_difference <- max(abs(
  as.matrix(summary_rows[, comparison_columns]) -
    as.matrix(authoritative_modelled[, comparison_columns])
))
if (largest_difference > 0.01) {
  stop(sprintf(
    "Rerun does not reproduce the packaged summaries (maximum difference %.4f).",
    largest_difference
  ))
}

saveRDS(
  list(
    model = "FishKelp v7 primary Albert-Chib ordered-probit model",
    target_species = target_species,
    routes = c("Trophic", "Habitat"),
    nchains = 4L,
    niter = 60000L,
    nburn = 15000L,
    max_split_rhat = max_rhat,
    draws = reach_draws
  ),
  file.path(output_dir, "blamey_bolton_fish_posterior_draws_v7.rds"),
  compress = "xz"
)

# Update the smaller draw file used by the focal Quarto report.
focal_species <- summary_rows$species[[1]]
focal_position <- match(focal_species, target_species)
focal_draws <- data.frame(
  chain = rep(seq_along(chains), each = n_saved_per_chain),
  iteration = rep(seq_len(n_saved_per_chain), times = length(chains)),
  Trophic = reach_draws[, focal_position, "Trophic"],
  Habitat = reach_draws[, focal_position, "Habitat"]
)
saveRDS(
  list(
    species = focal_species,
    model = "FishKelp v7 primary Albert-Chib ordered-probit model",
    niter = 60000L,
    nburn = 15000L,
    nchains = 4L,
    draws = focal_draws
  ),
  file.path(figure_data_dir, "strongest_fish_posterior_draws_v7.rds"),
  compress = "xz"
)

cat(sprintf(
  paste0(
    "Complete: %d posterior draws for each of %d fish on two routes. ",
    "Maximum split-Rhat %.3f; packaged summaries reproduced within %.4f.\n"
  ),
  n_draws,
  length(target_species),
  max_rhat,
  largest_difference
))
cat("Outputs -> ", output_dir, "\n", sep = "")
