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

grouped_class = scores_preprocessed_hierarchy.groupby("class")
grouping = "subclass"
for group in grouped_class.groups:
    grp = grouped_class.get_group(group).groupby(grouping).filter(lambda x: len(x) > 6)
    if len(grp) > 0 and group == "Benzene and substituted derivatives":
        fig = create_plot(grp, grouping, showlegend=False, hide_labels=True)
        fig.savefig(f"paper_plots/Fig6_benzene_subclasses.png", bbox_inches='tight')