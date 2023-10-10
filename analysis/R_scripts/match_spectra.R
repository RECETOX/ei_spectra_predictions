# Load required packages
library(MetaboAnnotation)
library(Spectra)
library(MsBackendMsp)
library(dplyr)
library(pheatmap)
library(ggplot2)


# Define file paths
reference_file <- '/home/wrojas/dev/RECETOX/ei_spectra_predictions/data/RECETOX_GC-EI_MS_20201028.msp'
simulated_file <- '/home/wrojas/dev/RECETOX/ei_spectra_predictions/data/matchms_filtering_default_filter_0.01_int_70_700_mz.msp'
output_file <- "/home/wrojas/dev/RECETOX/ei_spectra_predictions/analysis/R_scripts/output/matchSpectra_results_source_0.01_int.tsv"

# Read the reference and simulated spectra from MSP files
data_reference <- Spectra(reference_file, source = MsBackendMsp::MsBackendMsp())
data_simulated <- Spectra(simulated_file, source = MsBackendMsp::MsBackendMsp())


#################################
# matchSpectra
#################################
match_param <- MatchForwardReverseParam(requirePrecursor = FALSE,
    ppm = 5, 
    FUN = MsCoreUtils::ndotproduct, 
    THRESHFUN = function(x) which(x >= 0.0001), 
    THRESHFUN_REVERSE = function(x) which(x >= 0.0001) 
)

matched_spectra <- matchSpectra(data_simulated, data_reference, match_param)

print(matched_spectra)

matched_spectra_df <- spectraData(matched_spectra, columns = c("name", "target_name", "reverse_score", "score", "presence_ratio", "matched_peaks_count"))

matched_spectra_df <- as.data.frame(matched_spectra_df)

matched_spectra_df <- matched_spectra_df %>%
    filter_all(all_vars(!is.na(.)))

write.table(matched_spectra_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

#################################
# compareSpectra
#################################

# comp_spectra <- compareSpectra(data_simulated, data_reference, fun = "dotproduct", ppm = 5)

# output_file <- "/home/wrojas/dev/RECETOX/ei_spectra_predictions/analysis/R_scripts/output/comp_spectra.tsv"

# write.table(comp_spectra, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

# pheatmap_obj <- pheatmap(comp_spectra, cluster_rows = FALSE, cluster_cols = FALSE)

# ggsave("/home/wrojas/dev/RECETOX/ei_spectra_predictions/analysis/R_scripts/output/pheatmap.png", plot = pheatmap_obj, device = "png")




