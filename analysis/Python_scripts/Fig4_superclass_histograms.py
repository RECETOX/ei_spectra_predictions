import pandas as pd
import os
import numpy as np
import math
from matplotlib import pyplot as plt
from rdkit import Chem
import plotly.graph_objs as go
from plotly.subplots import make_subplots

from utils import *
from plotting import *

matchms_scores = load_matchms_scores()

matchms_scores_superclass = preprocess_data(normalize_df(matchms_scores.copy()), ["superclass"])
larger_superclasses = matchms_scores_superclass.groupby("superclass").filter(lambda x: len(x) > 2)
create_plot(larger_superclasses, "superclass", normalized_matches=True).savefig("paper_plots/Fig4a_superclasses_boxplot.png", bbox_inches='tight')
# plot name in the manuscript: "superclasses/20240207_boxplot_superclasses.png"

matches_normalized = matchms_scores['matches'] / matchms_scores['n_peaks_reference']
plt.clf()
plt.set_cmap('viridis')
plt.hist2d(matches_normalized * 100, matchms_scores['scores'] * 1000, bins=(5, 5), range=[[0, 100], [0, 1000]])
plt.colorbar()
plt.clim(0, 70)
plt.xlabel('ions matching reference (%)', fontsize=20)
plt.ylabel('scores', fontsize=20)
plt.tick_params(labelsize=13)
plt.gcf().set_size_inches(8, 6)
plt.savefig("paper_plots/Fig4a_superclasses_histogram.png", bbox_inches='tight')


matchms_scores_top5 = pd.read_csv("../data/output_matching/matchms/matchms_tol_0.0035_1%I_top5_with_0s_only_matching.tsv", sep="\t")
matchms_scores_top5.rename(columns={'CosineHungarian_0.0035_0.0_1.0_scores': 'scores'}, inplace=True)
matchms_scores_top5.rename(columns={'CosineHungarian_0.0035_0.0_1.0_matches': 'matches'}, inplace=True)
matchms_scores_top5 = append_classes(matchms_scores_top5, 'query')
matchms_scores_top5 = append_spectrum_metadata(matchms_scores_top5)

matchms_scores_superclass_top5 = preprocess_data(normalize_df(matchms_scores_top5.copy(), matches_norm_col=None), ["superclass"])
larger_superclasses_top5 = matchms_scores_superclass_top5.groupby("superclass").filter(lambda x: len(x) > 2)
create_plot(larger_superclasses_top5, "superclass", normalized_matches=False).savefig("paper_plots/Fig4b_superclasses_boxplot.png", bbox_inches='tight')
# plot name in the manuscript: "superclasses/20240223_boxplot_superclasses_top5.png"

plt.clf()
plt.set_cmap('viridis')
plt.hist2d(matchms_scores_top5['matches'], matchms_scores_top5['scores'] * 1000, bins=([0,1,2,3,4,5], 5))
plt.colorbar()
plt.clim(0, 70)

plt.xlabel('ion matches', fontsize=20)
plt.ylabel('scores', fontsize=20)
plt.tick_params(labelsize=13)
plt.gcf().set_size_inches(8, 6)
plt.savefig("paper_plots/Fig4b_superclasses_histogram.png", bbox_inches='tight')