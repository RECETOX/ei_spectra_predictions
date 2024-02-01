# Analysis

This folder contains the scripts used to generate the results presented in the paper.
This includes the boxplots for classes and chemical composition as well as some other plots that didn't make it into the paper, like some scatterplots and histograms.

It also contains the R scripts used to compute the NIST scores using the MetaboAnnotation package from the RforMassSpectrometry suite.

## Getting started

All packages needed to run the R scripts and Jupyter notebooks can be installed using the provided conda environment file.

```bash
conda env create -f analysis_environment.yaml
conda activate rmas
```

## Data
The data folder contains the data used to generate the plots.
This includes the generated as well as the predicted spectra after filtering, the structures and related metadata as well as the matching outputs and post-processed results.