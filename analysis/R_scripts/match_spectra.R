# Load required packages
library(Spectra)
library(MsBackendMsp)

# Define file paths
reference_file <- '/home/wrojas/dev/RECETOX/ei_spectra_predictions/data/RECETOX_GC-EI_MS_20201028.msp'
simulated_file <- '/home/wrojas/dev/RECETOX/ei_spectra_predictions/data/matchms_filtering_default_filter_0.1_int_70_700_mz.msp'

# Read the reference and simulated spectra from MSP files
data_reference <- Spectra(reference_file, source = MsBackendMsp::MsBackendMsp())
data_simulated <- Spectra(simulated_file, source = MsBackendMsp::MsBackendMsp())

# Print the variables (metadata) associated with the reference and simulated spectra
spectraVariables(data_reference)
spectraVariables(data_simulated)

# Extract the peaks data from the reference spectra
ref_peaks <- peaksData(data_reference)

# Compare the simulated spectra to itself with a ppm tolerance of 20
#comp_spectra_ppm <- compareSpectra(data_simulated, ppm = 20)

# Compare the simulated spectra to the reference spectra with a tolerance of 0.2
comp_spectra <- compareSpectra(data_simulated, data_reference, tolerance = 0.2)

# Match the simulated spectra to the reference spectra with a tolerance of 0.2
match_param <- CompareSpectraParam(tolerance = 0.2)
matched_spectra <- matchSpectra(data_simulated, data_reference, match_param)