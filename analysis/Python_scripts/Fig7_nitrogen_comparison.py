from utils import *
from plotting import boxplot_comparison

all_peaks_same = pd.read_csv("../data/output_matching/matchms/matchms_tol_0.0035_1%I_all_peaks_with_0s_only_matching.tsv", sep='\t')
all_peaks_same.rename(columns={'CosineHungarian_0.0035_0.0_1.0_scores': 'scores'}, inplace=True)
all_peaks_same.rename(columns={'CosineHungarian_0.0035_0.0_1.0_matches': 'matches'}, inplace=True)
all_peaks_same = append_spectrum_metadata(all_peaks_same)
merged_all_peaks_same = normalize_df(append_classes(all_peaks_same, 'query'))
mdf_comp = preprocess_data(merged_all_peaks_same, ["composition"])

baseline_cols= ['C,H', 'C,H,O', 'C,H,O,S', 'C,Cl,H,O', 'Br,C,H,O', 'C,Cl,H', 'C,Cl,H,O,S', 'C,Cl,F,H,O', 'C,H,O,P', 'C,H,O,P,S']
mdf_comp_baseline = mdf_comp.loc[mdf_comp['composition'].isin(baseline_cols)]
mdf_comp_baseline.sort_index(axis=1, inplace=True)

nitrogen_cols = ['C,H,N', 'C,H,N,O','C,H,N,O,S', 'C,Cl,H,N,O', 'Br,C,H,N,O', 'C,Cl,H,N', 'C,Cl,H,N,O,S', 'C,Cl,F,H,N,O','C,H,N,O,P', 'C,H,N,O,P,S']
mdf_comp_nitrogen = mdf_comp.loc[mdf_comp['composition'].isin(nitrogen_cols)]
mdf_comp_nitrogen.sort_index(axis=1, inplace=True)

boxplot_comparison(
    mdf_comp_baseline,
    baseline_cols,
    mdf_comp_nitrogen,
    nitrogen_cols,
    'scores',
    colors=['crimson', 'deepskyblue']
).savefig("paper_plots/Fig7_scores.png", bbox_inches='tight')

boxplot_comparison(
    mdf_comp_baseline,
    baseline_cols,
    mdf_comp_nitrogen,
    nitrogen_cols,
    'matches',
    colors=["darkgoldenrod", "yellow"],
).savefig("paper_plots/Fig7_matches.png", bbox_inches='tight')