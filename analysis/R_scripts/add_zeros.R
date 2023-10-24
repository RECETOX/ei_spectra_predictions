# Load required packages
library(MetaboAnnotation)
library(Spectra)
library(MsBackendMsp)
library(dplyr)

reference_file <- file.path('analysis/data/experimental/RECETOX_GC-EI_MS_20201028_norm_matchms.msp')
simulated_file <- file.path('analysis/data/filtered/matchms_filtering_default_filter_1%_all_peaks.msp')
matchms_match <- read.table(file.path('analysis/data/matchms_top5_comparison.tsv'), header = TRUE, sep = "\t", fill = TRUE)

# Check for NA values in the 'query' column
if (any(is.na(matchms_match$query))) {
  print("There are NA values in the 'query' column.")
}

# Check for NA values in the 'reference' column
if (any(is.na(matchms_match$reference))) {
  print("There are NA values in the 'reference' column.")
}

# Read the reference and simulated spectra from MSP files
data_reference <- Spectra(reference_file, source = MsBackendMsp::MsBackendMsp())
data_simulated <- Spectra(simulated_file, source = MsBackendMsp::MsBackendMsp())

sim_names <- data_simulated$name
ref_names <- data_reference$name

name_combination <- expand.grid(sim_names, ref_names)
zero_scores <- data.frame(query = character(), reference = character(), CosineHungarian_0.01_0.0_1.0_scores = numeric(), CosineHungarian_0.01_0.0_1.0_matches = numeric())

for (i in 1:nrow(name_combination)) {
  matches <- matchms_match[matchms_match$query == name_combination$Var1[i] & matchms_match$reference == name_combination$Var2[i],]
  if(nrow(matches) == 0) {
    new_row <- data.frame(name_combination$Var1[i], name_combination$Var2[i], 0, 0)
    names(new_row) <- c("query", "reference", "CosineHungarian_0.01_0.0_1.0_scores", "CosineHungarian_0.01_0.0_1.0_matches")
    zero_scores <- rbind(zero_scores, new_row)
  }
}

names(zero_scores) <- c("query", "reference", "CosineHungarian_0.01_0.0_1.0_scores", "CosineHungarian_0.01_0.0_1.0_matches")
matchms_match <- rbind(matchms_match, zero_scores)
write.table(matchms_match, file = "matchms_top5_comparison_with_zeros.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
