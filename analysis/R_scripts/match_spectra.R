# Load required packages
library(MetaboAnnotation)
library(Spectra)
library(MsBackendMsp)
library(dplyr)

match_spectra <- function(reference_file, simulated_file, output_file, ppm = 5) {
  # Read the reference and simulated spectra from MSP files
  data_reference <- Spectra(reference_file, source = MsBackendMsp::MsBackendMsp())
  data_simulated <- Spectra(simulated_file, source = MsBackendMsp::MsBackendMsp())


  # Define match parameters
  match_param <- MetaboAnnotation::MatchForwardReverseParam(
    requirePrecursor = FALSE,
    ppm = ppm, 
    FUN = MsCoreUtils::ndotproduct, 
    THRESHFUN = function(x) which(x >= 0.0) | TRUE, 
    THRESHFUN_REVERSE = function(x) which(x >= 0.0) | TRUE
  )

  # Perform matching
  matched_spectra <- MetaboAnnotation::matchSpectra(data_simulated, data_reference, match_param)

  # Convert matched spectra to data frame
  matched_spectra_df <- spectraData(matched_spectra, columns = c("name", "target_name", "reverse_score", "score", "presence_ratio", "matched_peaks_count"))
  matched_spectra_df <- as.data.frame(matched_spectra_df)
  write.table(matched_spectra_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

  # Write results to output file
  matched_spectra_df[is.na(matched_spectra_df[[3]]), 3] <- 0

  write.table(matched_spectra_df, file = paste(output_file,"filterNA",sep = '_'), sep = "\t", quote = FALSE, row.names = FALSE)
}

base_path <- 'analysis/data'

reference_file <- file.path(base_path, 'experimental/RECETOX_GC-EI_MS_20201028_matchms_filter_top5.msp')
simulated_file <- file.path(base_path, 'filtered/1%/Galaxy29-[matchms_filtering_on_data_27]_1%_tp5.msp')
match_spectra(
  reference_file,
  simulated_file,
  file.path(base_path, 'outputs/matchspectra_R/metaboannotation_top5.tsv')
)

# match_spectra(
#   file.path(base_path, "experimental/RECETOX_GC-EI_MS_20201028_norm_matchms.msp"),
#   file.path(base_path, 'filtered/matchms_filtering_default_filter_1%_all_peaks.msp'),
#   file.path(base_path, 'outputs/matchspectra_R/match_1%/metaboannotation_all_peaks_zeros_temp.tsv')
# )

# reference_file <- file.path(base_path, "experimental/problematic_spectra.msp")
# simulated_file <- file.path(base_path, 'filtered/problematic_predicted.msp')
# match_spectra(
#   reference_file,
#   query_file,
#   file.path(base_path, 'outputs/matchspectra_R/problematic_scores.tsv')
# )
# match_spectra(
#   file.path(base_path, "experimental/RECETOX_GC-EI_MS_20201028.msp"),
#   file.path(base_path, 'filtered/matchms_filtering_default_filter_1%_all_peaks.msp'),
#   file.path(base_path, 'outputs/matchspectra_R/match_1%/matchSpectra_results_source_1%_all_peaks.tsv')
# )