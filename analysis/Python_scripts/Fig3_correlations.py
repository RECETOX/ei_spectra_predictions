import pandas as pd
from matplotlib import pyplot as plt

from utils import *
from plotting import *

matchms_scores = load_matchms_scores()
_, _, label_fontsize, tick_fontsize, text_width, _ = init()

df = normalize_df(matchms_scores, matches_norm_col=None)
del df['peak_comments']

matches_col = 'matches'
scores_col = 'scores'

df['matches_norm_query'] = df[matches_col] / df['n_peaks_query']
df['matches_norm_reference'] = df[matches_col] / df['n_peaks_reference']

properties = [
    'scores',
    'matches',
    'matches_norm_query',
    'matches_norm_reference',
    'molecular_flexibility',
    'rotatable_bonds',
    'stereo_centers',
    'molecular_complexity',
    'n_atoms',
    'precursor_mz',
    'electronegative_atoms',
    'aromatic_nitrogens',
    'amines',
    'amides',
]

# Assuming `df` is your DataFrame
corr = df[properties].corr().round(2)

plt.figure(figsize=(24, 20))
cax = sns.heatmap(
    corr,
    annot=True,
    cmap='coolwarm',
    center=0,
    vmin=-1,
    vmax=1,
    annot_kws={"size": label_fontsize}
)
# plt.title('Pearson Correlations')
plt.tick_params(axis='both', which='major', labelsize=tick_fontsize)
# Get the colorbar from the HeatMap and set the fontsize for its tick labels
cbar = cax.collections[0].colorbar
cbar.ax.tick_params(labelsize=tick_fontsize)

plt.savefig("paper_plots/Fig3_correlations.png", bbox_inches='tight')
# plot name in the manuscript:
# "correlations/20240517_heatmap_properties_correlations.png"