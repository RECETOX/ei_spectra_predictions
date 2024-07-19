# %%
from matchms.importing import load_from_msp
from matchms.plotting import plot_spectra_mirror
import pandas as pd
from matplotlib import pyplot as plt
from pathlib import Path

# %%
scores = pd.read_csv("../data/output_matching/matchms/matchms_tol_0.0035_1%I_all_peaks_with_0s_only_matching.tsv", sep='\t')

# %%
query_spectra = list(load_from_msp("../data/filtered/simulated_matchms_filter_1%I_all_peaks.msp"))
query_lookup = { x.get("compound_name"): x for x in query_spectra }

# %%
reference_spectra = list(load_from_msp("../data/experimental/RECETOX_GC-EI_MS_20201028.msp"))
reference_lookup = { x.get("compound_name"): x for x in reference_spectra }

# %%
if not Path.exists(Path("mirrorplots")):
    Path.mkdir("mirrorplots")

for idx, row in scores.iterrows():
    key_query:str  = row["query"]
    key_ref:str  = row["reference"]

    query_spectrum = query_lookup.get(key_query)
    reference_spectrum = reference_lookup.get(key_ref)

    if query_spectrum and reference_spectrum:
        fig, ax = plt.subplots(figsize=(12, 4))
        plot_spectra_mirror(query_spectrum, reference_spectrum, ax=ax, color_top="red", color_bottom="blue")

        label = key_query.replace("/","-")

        fig.savefig(f"mirrorplots/{idx}_{label}.png", bbox_inches='tight')
        plt.close()
        plt.clf()


