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

##############################################
# Simulated spectra must be in MS-LIMA format.
##############################################

# matching for all peaks
reference_file_all_peaks <- file.path('analysis/data/experimental/RECETOX_GC-EI_MS_20201028.msp')
simulated_file_all_peaks <- file.path('analysis/data/filtered/simulated_matchms_filter_1%I_all_peaks.msp')

match_spectra(
  reference_file_all_peaks,
  simulated_file_all_peaks,
  file.path('analysis/data/output_matching/metaboannotation/metaboannotation_NA_to_zeros_all_peaks.tsv')
)

# matching for top5 peaks. query and reference were filtered to top5 peaks.
reference_file_top5 <- file.path('analysis/data/experimental/RECETOX_GC-EI_MS_20201028_matchms_filter_top5.msp')
simulated_file_top5 <- file.path('analysis/data/filtered/simulated_matchms_filter_1%I_top5.msp')

match_spectra(
  reference_file_top5,
  simulated_file_top5,
  file.path('analysis/data/output_matching/metaboannotation/metaboannotation_NA_to_zeros_top5.tsv')
)