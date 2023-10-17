# Load required packages
library(MetaboAnnotation)
library(Spectra)
library(MsBackendMsp)
library(dplyr)

match_spectra <- function(reference_file, simulated_file, output_file, ppm = 5, requirePrecursor = FALSE, FUN = MsCoreUtils::ndotproduct, THRESHFUN = function(x) which(x >= 0.0001), THRESHFUN_REVERSE = function(x) which(x >= 0.0001)) {
  # Read the reference and simulated spectra from MSP files
  data_reference <- Spectra(reference_file, source = MsBackendMsp::MsBackendMsp())
  data_simulated <- Spectra(simulated_file, source = MsBackendMsp::MsBackendMsp())

  # Define match parameters
  match_param <- MatchForwardReverseParam(requirePrecursor = requirePrecursor,
                                          ppm = ppm, 
                                          FUN = FUN, 
                                          THRESHFUN = THRESHFUN, 
                                          THRESHFUN_REVERSE = THRESHFUN_REVERSE)

  # Perform matching
  matched_spectra <- matchSpectra(data_simulated, data_reference, match_param)

  # Convert matched spectra to data frame
  matched_spectra_df <- spectraData(matched_spectra, columns = c("name", "target_name", "reverse_score", "score", "presence_ratio", "matched_peaks_count"))
  matched_spectra_df <- as.data.frame(matched_spectra_df)
  matched_spectra_df <- matched_spectra_df %>%
    filter_all(all_vars(!is.na(.)))

  # Write results to output file
  write.table(matched_spectra_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

reference_file <- '/home/wrojas/dev/RECETOX/ei_spectra_predictions/data/RECETOX_GC-EI_MS_20201028_norm_matchms.msp'
query_file <- '/home/wrojas/dev/RECETOX/ei_spectra_predictions/analysis/data/filtered/matchms_filtering_default_filter_1%_all_peaks.msp'
output_file <- "/home/wrojas/dev/RECETOX/ei_spectra_predictions/analysis/data/matching/matchspectra_R/match_1%/matchSpectra_results_source_1%_all_peaks.tsv"

match_spectra(reference_file, query_file, output_file)