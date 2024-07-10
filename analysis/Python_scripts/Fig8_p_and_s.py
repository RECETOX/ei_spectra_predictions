from utils import *
from plotting import create_plot

matchms_scores = load_matchms_scores()
merged_all_peaks_same = normalize_df(matchms_scores)
mdf_comp = preprocess_data(merged_all_peaks_same, ["composition"])

mdf_comp_ps = mdf_comp[mdf_comp['composition'].str.contains('S|P')]
mdf_comp_ps = mdf_comp_ps[mdf_comp_ps['composition'] != 'C,F,H,N,Si']
mdf_comp_ps = mdf_comp_ps.groupby('composition').filter(lambda x: len(x) > 2)
create_plot(mdf_comp_ps, "composition").savefig("paper_plots/Fig8_P_and_S.png", bbox_inches='tight')