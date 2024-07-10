import pandas as pd
from matchms.logging_functions import set_matchms_logger_level

from utils import append_classes, load_spectra_metadata, normalize_df
from plotting import scatterplot_matplotlib

set_matchms_logger_level('ERROR')

matchms_scores = pd.read_csv("../data/output_matching/matchms/matchms_tol_0.0035_1%I_all_peaks_with_0s_only_matching.tsv", sep="\t")
matchms_scores.rename(columns={'CosineHungarian_0.0035_0.0_1.0_scores': 'scores'}, inplace=True)
matchms_scores.rename(columns={'CosineHungarian_0.0035_0.0_1.0_matches': 'matches'}, inplace=True)

_ , spectra_metadata, _ = load_spectra_metadata("../data/filtered/simulated_matchms_filter_1%I_all_peaks.msp", 'query')
_ , reference_spectra_metadata, _ = load_spectra_metadata("../data/experimental/RECETOX_GC-EI_MS_20201028.msp", 'reference')

merged = matchms_scores.merge(spectra_metadata, on="query", how="inner")
merged.rename(columns={'num_peaks': 'n_peaks_query'}, inplace=True)

merged = merged.merge(reference_spectra_metadata, on="reference", how="inner")
merged.rename(columns={'num_peaks': 'n_peaks_reference'}, inplace=True)

numeric_columns = ['matches', 'n_peaks_query', 'n_peaks_reference']
merged[numeric_columns] = merged[numeric_columns].apply(pd.to_numeric, errors='coerce')

merged['FractionQuery'] = merged['matches'] / merged['n_peaks_query']
merged['FractionReference'] = merged['matches'] / merged['n_peaks_reference']

merged = append_classes(merged, "query")

# Create a scatter plot
scatterplot_matplotlib(normalize_df(merged, matches_norm_col=None)).savefig("paper_plots/Fig2_scatterplot.png", bbox_inches='tight')
# plot name in the manuscript:
# "20240517_scatterplot.png"