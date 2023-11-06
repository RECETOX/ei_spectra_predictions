# Load required packages
library(MetaboAnnotation)
library(Spectra)
library(MsBackendMsp)
library(dplyr)


add_zeros <- function(out_file, match_file, reference_file, simulated_file, query, reference, scores, matches) {
  # if (any(is.na(match_file$query))) {
  #   print("There are NA values in the 'query' column.")
  # }

  # # Check for NA values in the 'reference' column
  # if (any(is.na(match_file$reference))) {
  #   print("There are NA values in the 'reference' column.")
  # }

  # Read the reference and simulated spectra from MSP files
  data_reference <- Spectra(reference_file, source = MsBackendMsp::MsBackendMsp())
  data_simulated <- Spectra(simulated_file, source = MsBackendMsp::MsBackendMsp())

  sim_names <- data_simulated$name
  ref_names <- data_reference$name
  #match_file[is.na(match_file[[3]]), 3] <- 0

  name_combination <- expand.grid(sim_names, ref_names)
  zero_scores <- data.frame(matrix(ncol = ncol(match_file), nrow = 0))
  names(zero_scores) <- names(match_file)

  for (i in 1:nrow(name_combination)) {
    matches <- match_file[match_file[query] == as.character(name_combination$Var1[i]) & match_file[reference] == as.character(name_combination$Var2[i]),]
    if(nrow(matches) == 0) {
       new_row <- data.frame(name_combination$Var1[i], name_combination$Var2[i], 0, 0)
       zero_scores <- rbind(zero_scores, new_row)
    }
  }
  names(zero_scores) <- names(match_file)

  match_file <- rbind(match_file, zero_scores)
  write.table(match_file, file = out_file)
}


reference_file <- file.path('analysis/data/experimental/RECETOX_GC-EI_MS_20201028_norm_matchms.msp')
simulated_file <- file.path('analysis/data/filtered/matchms_filtering_default_filter_1%_all_peaks.msp')

# for MatchMS top 5
matchms_match <- read.table(file.path('analysis/data/matchms_top5_comparison.tsv'), header = TRUE, sep = "\t", fill = TRUE)
out_file_matchms = file.path('analysis/data/matchms_top5_comparison_with_zeros.tsv')
add_zeros(out_file_matchms, matchms_match, reference_file, simulated_file, 
 "query", "reference", "CosineHungarian_0.01_0.0_1.0_scores", "CosineHungarian_0.01_0.0_1.0_matches")

# for MatchMS all peaks
matchms_match <- read.table(file.path('analysis/data/matchms_all_peaks.tsv'), header = TRUE, sep = "\t", fill = TRUE)
out_file_matchms = file.path('analysis/data/matchms_all_peaks_with_zeros.tsv')
add_zeros(out_file_matchms, matchms_match, reference_file, simulated_file, 
 "query", "reference", "CosineHungarian_0.01_0.0_1.0_scores", "CosineHungarian_0.01_0.0_1.0_matches")

# # for MatchSpectra
# mspectra_match <- read.table(file.path('analysis/data/outputs/matchspectra_R/metaboannotation_top5.tsv'), header = TRUE, sep = "\t", fill = TRUE)
# out_file_mspectra = file.path('analysis/R_scripts/metaboannotation_top5_with_zeros.tsv')
# add_zeros(out_file_mspectra, mspectra_match, reference_file, simulated_file,
#  "name", "target_name", "score", "matched_peaks_count")