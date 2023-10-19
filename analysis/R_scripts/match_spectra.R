# Load required packages
library(MetaboAnnotation)
library(Spectra)
library(MsBackendMsp)
library(dplyr)

match_spectra <- function(reference_file, simulated_file, output_file, ppm = 5, 
                          requirePrecursor = FALSE, FUN = MsCoreUtils::ndotproduct,
                          THRESHFUN = function(x) which(x >= 0.0), THRESHFUN_REVERSE = function(x) which(x >= 0.0)) {
  # Read the reference and simulated spectra from MSP files
  data_reference <- Spectra(reference_file, source = MsBackendMsp::MsBackendMsp())
  data_simulated <- Spectra(simulated_file, source = MsBackendMsp::MsBackendMsp())


  # Define match parameters
  match_param <- MetaboAnnotation::MatchForwardReverseParam(requirePrecursor = requirePrecursor,
                                          ppm = ppm, 
                                          FUN = FUN, 
                                          THRESHFUN = THRESHFUN, 
                                          THRESHFUN_REVERSE = THRESHFUN_REVERSE)

  # Perform matching
  matched_spectra <- MetaboAnnotation::matchSpectra(data_simulated, data_reference, match_param)

  # Convert matched spectra to data frame
  matched_spectra_df <- spectraData(matched_spectra, columns = c("name", "target_name", "reverse_score", "score", "presence_ratio", "matched_peaks_count"))
  matched_spectra_df <- as.data.frame(matched_spectra_df)
 # matched_spectra_df <- matched_spectra_df # %>%    filter_all(all_vars(!is.na(.)))

  # matched_spectra_filtered <- subset(matched_spectra_df, "name" == "target_name")
  # missing_names <- setdiff(data_simulated$name, matched_spectra_filtered$name)
  # print(missing_names)
  # file_path <- "missing_names_nofilter_NA.txt"
  # write.table(missing_names, file_path, col.names = FALSE, row.names = FALSE, quote = FALSE)
 
  # Write results to output file
  write.table(matched_spectra_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)
}

base_path <- '/home/wrojas/dev/RECETOX/ei_spectra_predictions/analysis/data'

match_spectra(
  file.path(base_path, "experimental/RECETOX_GC-EI_MS_20201028_norm_matchms.msp"),
  file.path(base_path, 'filtered/matchms_filtering_default_filter_1%_all_peaks.msp'),
  file.path(base_path, 'outputs/matchspectra_R/match_1%/metaboannotation_all_peaks_zeros_temp.tsv')
)

# match_spectra(
#   file.path(base_path, "experimental/RECETOX_GC-EI_MS_20201028.msp"),
#   file.path(base_path, 'filtered/matchms_filtering_default_filter_1%_all_peaks.msp'),
#   file.path(base_path, 'outputs/matchspectra_R/match_1%/matchSpectra_results_source_1%_all_peaks.tsv')
# )