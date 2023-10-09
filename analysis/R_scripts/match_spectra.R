# Load required packages
library(MetaboAnnotation)
library(Spectra)
library(MsBackendMsp)
library(dplyr)
library(pheatmap)
library(ggplot2)

# Define file paths
reference_file <- '/home/wrojas/dev/RECETOX/ei_spectra_predictions/data/RECETOX_GC-EI_MS_20201028.msp'
simulated_file <- '/home/wrojas/dev/RECETOX/ei_spectra_predictions/data/matchms_filtering_default_filter_0.1_int_70_700_mz.msp'
output_file <- "/home/wrojas/dev/RECETOX/ei_spectra_predictions/analysis/R_scripts/output/matched_spectra_df.tsv"

# Read the reference and simulated spectra from MSP files
data_reference <- Spectra(reference_file, source = MsBackendMsp::MsBackendMsp())
data_simulated <- Spectra(simulated_file, source = MsBackendMsp::MsBackendMsp())


#################################
# matchSpectra
#################################
# Match the simulated spectra to the reference spectra using matchSpectra function
match_param <- MatchForwardReverseParam(
    ppm = 10, # Set the mass tolerance to 10 ppm 
    FUN = MsCoreUtils::ndotproduct, # Use the dot product as the similarity function
    THRESHFUN = function(x) which(x >= 0.01), # Set the threshold function for forward matching
    THRESHFUN_REVERSE = function(x) which(x >= 0.01) # Set the threshold function for reverse matching
)
matched_spectra <- matchSpectra(data_simulated, data_reference, match_param)

# Print the matched spectra object
print(matched_spectra)

# Convert the matched spectra to a data frame
matched_spectra_df <- spectraData(matched_spectra, columns = c("name", "target_name", "reverse_score", "score", "presence_ratio", "matched_peaks_count"))

# Convert the matched spectra data frame to a standard data frame
matched_spectra_df <- as.data.frame(matched_spectra_df)

# Replace NA values with 0.0 in the matched spectra data frame using dplyr
matched_spectra_df <- matched_spectra_df %>%
    mutate_all(~ifelse(is.na(.), 0.0, .))

# Write the matched spectra data frame to a TSV file
write.table(matched_spectra_df, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

#################################
# compareSpectra
#################################

# Match the spectra using compareSpectra function 
comp_spectra <- compareSpectra(data_simulated, data_reference, fun = "dotproduct", ppm = 5)

# Define the output file path for the comparison results
output_file <- "/home/wrojas/dev/RECETOX/ei_spectra_predictions/analysis/R_scripts/output/comp_spectra.tsv"

# Write the comparison results to the output TSV file
write.table(comp_spectra, file = output_file, sep = "\t", quote = FALSE, row.names = FALSE)

# Create a pheatmap object from the comparison results
pheatmap_obj <- pheatmap(comp_spectra, cluster_rows = FALSE, cluster_cols = FALSE)

# Save the pheatmap object to a PNG file
ggsave("/home/wrojas/dev/RECETOX/ei_spectra_predictions/analysis/R_scripts/output/pheatmap.png", plot = pheatmap_obj, device = "png")