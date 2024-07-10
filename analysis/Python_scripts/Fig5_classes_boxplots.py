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
merged = normalize_df(matchms_scores.copy())

scores_preprocessed_hierarchy = preprocess_data(merged, ["superclass", "class", "subclass"])
grouped_superclass = scores_preprocessed_hierarchy.groupby("superclass")
grouping = "class"

for group in grouped_superclass.groups:
    grp = grouped_superclass.get_group(group).groupby(grouping).filter(lambda x: len(x) > 2)
    if len(grp) > 0:
        fig = create_plot(grp, grouping, showlegend=False, hide_labels=True)
        fig.savefig(f"paper_plots/Fig5_{group}.png", bbox_inches='tight')
# plot name in the manuscript in that order:
# "classes/20240207_boxplot_benzenoids.png"
# "classes/20240207_boxplot_lipids.png"
# "classes/20240207_boxplot_organic_acids.png"
# "classes/20240207_boxplot_organooxygen.png"
# "classes/20240207_boxplot_organohalogen.png"
# "classes/20240207_boxplot_organoheterocyclic.png"
# "classes/20240207_boxplot_phenylpropanoids.png"

