library(MetaboAnnotation)
library(Spectra)
library(MsBackendMsp)
# read the MSP file  into a data frame
data_reference <- system.file("extdata", "/home/wrojas/dev/RECETOX/ei_spectra_predictions/data/RECETOX_GC-EI_MS_20201028.msp", package = "MsBackendMsp")
dr <- Spectra(data_reference, source = MsBackendMsp())
spectraVariables(dr)

# data_simulated <- system.file("extdata", "/home/wrojas/dev/RECETOX/ei_spectra_predictions/data/all_results.msp", package = "MsBackendMsp")
# ds <- Spectra(data_simulated, source = MsBackendMsp())

