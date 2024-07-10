import numpy as np
import math
from matplotlib import pyplot as plt
from rdkit import Chem
import plotly.graph_objs as go
from plotly.subplots import make_subplots

from utils import *
from plotting import *


matchms_scores = pd.read_csv("../data/output_matching/matchms/matchms_tol_0.0035_1%I_all_peaks_with_0s_only_matching.tsv", sep="\t")
matchms_scores.rename(columns={'CosineHungarian_0.0035_0.0_1.0_scores': 'scores'}, inplace=True)
matchms_scores.rename(columns={'CosineHungarian_0.0035_0.0_1.0_matches': 'matches'}, inplace=True)
matchms_scores = append_classes(matchms_scores, 'query')
matchms_scores = append_spectrum_metadata(matchms_scores)
merged = normalize_df(matchms_scores.copy())

scores_preprocessed_hierarchy = preprocess_data(merged, ["superclass", "class", "subclass"])

grouped_class = scores_preprocessed_hierarchy.groupby("class")
grouping = "subclass"
for group in grouped_class.groups:
    grp = grouped_class.get_group(group).groupby(grouping).filter(lambda x: len(x) > 6)
    if len(grp) > 0 and group == "Benzene and substituted derivatives":
        fig = create_plot(grp, grouping, showlegend=False, hide_labels=True)
        fig.savefig(f"paper_plots/Fig6_benzene_subclasses.png", bbox_inches='tight')