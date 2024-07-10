from utils import *
from plotting import create_plot

all_peaks_same = pd.read_csv("../data/output_matching/matchms/matchms_tol_0.0035_1%I_all_peaks_with_0s_only_matching.tsv", sep='\t')
all_peaks_same.rename(columns={'CosineHungarian_0.0035_0.0_1.0_scores': 'scores'}, inplace=True)
all_peaks_same.rename(columns={'CosineHungarian_0.0035_0.0_1.0_matches': 'matches'}, inplace=True)
all_peaks_same = append_spectrum_metadata(all_peaks_same)
merged_all_peaks_same = normalize_df(append_classes(all_peaks_same, 'query'))
mdf_comp = preprocess_data(merged_all_peaks_same, ["composition"])

mdf_comp_ps = mdf_comp[mdf_comp['composition'].str.contains('S|P')]
mdf_comp_ps = mdf_comp_ps[mdf_comp_ps['composition'] != 'C,F,H,N,Si']
mdf_comp_ps = mdf_comp_ps.groupby('composition').filter(lambda x: len(x) > 2)
create_plot(mdf_comp_ps, "composition").savefig("paper_plots/Fig8_P_and_S.png", bbox_inches='tight')