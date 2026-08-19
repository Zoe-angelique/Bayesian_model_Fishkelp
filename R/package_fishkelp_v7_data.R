#!/usr/bin/env Rscript

# Package the FishKelp v7 data needed for the fish listed in Table 2 of
# Blamey & Bolton (2018). Run from the Bayesian_model_Fishkelp repository root:
#
#   Rscript R/package_fishkelp_v7_data.R /path/to/FishKelp

fishkelp_root <- normalizePath(
  commandArgs(trailingOnly = TRUE)[[1]],
  mustWork = TRUE
)

repo_root <- normalizePath(".", mustWork = TRUE)
output_dir <- file.path(repo_root, "data", "fishkelp_v7")
model_output_dir <- file.path(output_dir, "model")
dir.create(model_output_dir, recursive = TRUE, showWarnings = FALSE)

kelp_project <- file.path(fishkelp_root, "projects", "kelp-dependence")
source_paths <- c(
  coded_records = file.path(
    kelp_project, "data", "rfbf_v7", "coded_responses_for_model_v7.csv"
  ),
  analysis_records = file.path(
    kelp_project, "data", "rfbf_v7", "analysis_dataset_v7.csv"
  ),
  species_frame = file.path(
    kelp_project, "data", "rfbf_v7", "species_frame_v7.csv"
  ),
  species_traits = file.path(
    kelp_project, "metadata", "species_traits_v7.csv"
  ),
  posterior_estimates = file.path(
    kelp_project, "data", "model_v7", "tdr_hdr_estimates_v7.csv"
  ),
  model_data = file.path(
    kelp_project, "data", "model_v7", "tdr_model_data_v7.rds"
  ),
  convergence = file.path(
    kelp_project, "data", "model_v7", "tdr_convergence_v7.csv"
  ),
  sampler = file.path(
    kelp_project, "code", "model_v7", "tdr_sampler_v7.R"
  )
)

# Table 2 as published. accepted_name records the name used by FishKelp v7.
paper_fish <- data.frame(
  table_order = 1:30,
  common_name_published = c(
    "Baardman", "Black musselcracker", "Blacktail", "Bronze bream",
    "Cape knifejaw", "Cape stumpnose", "Carpenter",
    "Common smooth-hound shark", "Dageraad", "Dorado", "Elf",
    "Englishman", "Garrick", "Geelbek", "Hottentot", "Janbruin",
    "Panga", "Roman", "Red steenbras", "Red stumpnose", "Santer",
    "Scotsman", "Snoek", "Soupfin shark", "Pyjama shark",
    "White musselcracker", "White stumpnose", "Yellow-belly rockcod",
    "Yellowtail", "Zebra fish"
  ),
  scientific_name_published = c(
    "Umbrina sp.", "Cymatoceps nasutus", "Diplodus sargus capensis",
    "Pachymetopon grande", "Oplegnathus conwayi", "Rhabdosargus holubi",
    "Argyrozona argyrozona", "Mustelus mustelus",
    "Chrysoblephus cristiceps", "Coryphaena hippurus",
    "Pomatomus saltatrix", "Chrysoblephus anglicus", "Lichia amia",
    "Atractoscion aequidens", "Pachymetopon blochii",
    "Gymnocrotaphus curvidens", "Pterogymnus laniarius",
    "Chrysoblephus laticeps", "Petrus rupestris",
    "Chrysoblephus gibbiceps", "Cheimerius nufar",
    "Polysteganus praeorbitalis", "Thyrsites atun", "Galeorhinus galeus",
    "Poroderma africanum", "Sparodon durbanensis",
    "Rhabdosargus globiceps", "Epinephelus marginatus", "Seriola lalandi",
    "Diplodus cervinus hottentotus"
  ),
  accepted_name = c(
    "Umbrina robinsoni", "Cymatoceps nasutus", "Diplodus capensis",
    "Pachymetopon grande", "Oplegnathus conwayi", "Rhabdosargus holubi",
    "Argyrozona argyrozona", "Mustelus mustelus",
    "Chrysoblephus cristiceps", "Coryphaena hippurus",
    "Pomatomus saltatrix", "Chrysoblephus anglicus", "Lichia amia",
    "Atractoscion aequidens", "Pachymetopon blochii",
    "Gymnocrotaphus curvidens", "Pterogymnus laniarius",
    "Chrysoblephus laticeps", "Petrus rupestris",
    "Chrysoblephus gibbiceps", "Cheimerius nufar",
    "Polysteganus praeorbitalis", "Thyrsites atun", "Galeorhinus galeus",
    "Poroderma africanum", "Sparodon durbanensis",
    "Rhabdosargus globiceps", "Epinephelus marginatus", "Seriola lalandi",
    "Diplodus hottentotus"
  ),
  fishery_published = c(
    "Recreational", "Recreational", "Recreational", "Recreational",
    "Recreational", "Recreational", "Commercial & recreational",
    "Commercial", "Commercial & recreational", "Commercial",
    "Commercial & recreational", "Commercial & recreational",
    "Recreational", "Commercial & recreational",
    "Commercial & recreational", "Recreational", "Commercial",
    "Commercial & recreational", "None (previously recreational)",
    "Commercial & recreational", "Commercial & recreational",
    "Commercial & recreational", "Commercial & recreational", "Commercial",
    "Recreational", "Recreational", "Commercial & recreational",
    "Commercial & recreational", "Commercial & recreational", "Recreational"
  ),
  exploitation_status_published = c(
    "Susceptible to overfishing",
    "Unknown, but vulnerable to fishing pressure",
    "Susceptible to overfishing", "IUCN near threatened", "Unknown",
    "IUCN least concern", "Optimally exploited (previously overfished)",
    "IUCN vulnerable", "Overexploited", "Unknown", "Unknown", "Collapsed",
    "Unknown", "Overexploited", "Sustainable", "Unknown", "Unknown",
    "Previously depleted", "Overexploited", "Overexploited", "Unknown",
    "Overexploited", "Optimally/fully exploited", "Fully/overexploited",
    "IUCN near threatened", "IUCN vulnerable", "IUCN vulnerable",
    "Unknown/IUCN endangered", "Sustainable", "Unknown"
  ),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

read_snapshot <- function(path) {
  read.csv(
    path,
    stringsAsFactors = FALSE,
    check.names = FALSE,
    na.strings = c("NA", "")
  )
}

coded_records <- read_snapshot(source_paths[["coded_records"]])
analysis_records <- read_snapshot(source_paths[["analysis_records"]])
species_frame <- read_snapshot(source_paths[["species_frame"]])
species_traits <- read_snapshot(source_paths[["species_traits"]])
posterior_estimates <- read_snapshot(source_paths[["posterior_estimates"]])
convergence <- read_snapshot(source_paths[["convergence"]])

paper_names <- paper_fish$accepted_name
paper_fish$in_species_frame_v7 <- paper_names %in% species_frame$species
scope_values <- species_frame$rfbf_scope[match(paper_names, species_frame$species)]
paper_fish$rfbf_scope_v7 <- NA
known_scope <- !is.na(scope_values)
paper_fish$rfbf_scope_v7[known_scope] <-
  tolower(as.character(scope_values[known_scope])) == "true"
paper_fish$modelled_v7 <- paper_names %in% posterior_estimates$species
paper_fish$n_model_records_v7 <- as.integer(table(factor(
  coded_records$species,
  levels = paper_names
)))
paper_fish$n_analysis_records_v7 <- as.integer(table(factor(
  analysis_records$species,
  levels = paper_names
)))

paper_fish$model_status_v7 <- "no_model_ready_records"
paper_fish$model_status_v7[
  is.na(paper_fish$rfbf_scope_v7) | !paper_fish$rfbf_scope_v7
] <- "outside_rfbf_scope"
paper_fish$model_status_v7[!paper_fish$in_species_frame_v7] <-
  "not_in_rfbf_species_frame"
paper_fish$model_status_v7[paper_fish$modelled_v7] <- "modelled"

paper_fish$source_table_note <- paste(
  "Table 2 pools kelp forest and temperate reef ecosystems;",
  "listing alone is not species-specific evidence of a kelp relationship."
)

add_paper_keys <- function(frame) {
  index <- match(frame$species, paper_fish$accepted_name)
  cbind(
    paper_table_order = paper_fish$table_order[index],
    paper_common_name = paper_fish$common_name_published[index],
    frame
  )
}

model_records_subset <- coded_records[coded_records$species %in% paper_names,]
model_records_subset <- add_paper_keys(model_records_subset)
model_records_subset <- model_records_subset[
  order(
    model_records_subset$paper_table_order,
    model_records_subset$route,
    model_records_subset$source_id
  ),
]

analysis_records_subset <- analysis_records[
  analysis_records$species %in% paper_names,
]
analysis_records_subset <- add_paper_keys(analysis_records_subset)
analysis_records_subset <- analysis_records_subset[
  order(
    analysis_records_subset$paper_table_order,
    analysis_records_subset$route,
    analysis_records_subset$source_id
  ),
]

left_join_paper <- function(frame, key = "species") {
  index <- match(paper_fish$accepted_name, frame[[key]])
  extra <- frame[index, setdiff(names(frame), key), drop = FALSE]
  rownames(extra) <- NULL
  cbind(paper_fish, extra)
}

posterior_subset <- left_join_paper(posterior_estimates)
species_frame_subset <- left_join_paper(species_frame)
traits_subset <- left_join_paper(species_traits)

write_snapshot <- function(frame, filename) {
  write.csv(
    frame,
    file.path(output_dir, filename),
    row.names = FALSE,
    na = "",
    fileEncoding = "UTF-8"
  )
}

write_snapshot(paper_fish, "blamey_bolton_2018_fish.csv")
write_snapshot(model_records_subset, "fishkelp_model_records_v7.csv")
write_snapshot(analysis_records_subset, "fishkelp_analysis_records_v7.csv")
write_snapshot(posterior_subset, "fishkelp_posterior_estimates_v7.csv")
write_snapshot(species_frame_subset, "fishkelp_species_frame_v7.csv")
write_snapshot(traits_subset, "fishkelp_species_traits_v7.csv")
write_snapshot(convergence, file.path("model", "tdr_convergence_v7.csv"))

# The full model-ready object is intentionally retained. The hierarchical model
# must be fitted to the complete FishKelp data and filtered to the paper species
# only after fitting.
saveRDS(
  readRDS(source_paths[["model_data"]]),
  file.path(model_output_dir, "tdr_model_data_v7.rds"),
  compress = "xz"
)

git_commit <- paste(tryCatch(
  system2(
    "git",
    c("-C", shQuote(fishkelp_root), "rev-parse", "HEAD"),
    stdout = TRUE,
    stderr = FALSE
  ),
  error = function(e) NA_character_
), collapse = "")

provenance <- data.frame(
  source_item = names(source_paths),
  fishkelp_relative_path = sub(
    paste0("^", fishkelp_root, "/?"),
    "",
    unname(source_paths)
  ),
  source_md5 = unname(tools::md5sum(source_paths)),
  fishkelp_git_commit = git_commit[[1]],
  packaged_on = as.character(Sys.Date()),
  stringsAsFactors = FALSE
)
write_snapshot(provenance, "SOURCE_PROVENANCE.csv")

cat(sprintf(
  paste0(
    "Packaged %d Blamey & Bolton fish: %d modelled, ",
    "%d model records and %d detailed analysis records.\n"
  ),
  nrow(paper_fish),
  sum(paper_fish$modelled_v7),
  nrow(model_records_subset),
  nrow(analysis_records_subset)
))
cat("Output -> ", output_dir, "\n", sep = "")
