import pandas as pd

from utils import *
from plotting import make_simple_boxplot

wang_db = pd.read_csv("../data/reference/13321_2020_470_MOESM1_ESM.csv")
wang_db_tms = pd.read_csv("../data/reference/wang2022_TMS_annotated.tsv", sep="\t", decimal=",")
df = normalize_df(load_matchms_scores())

wang_db["study"] = "Wang et al.$^{27}$"
wang_db_tms["study"] = "Wang et al.$^{28}$"
df["study"] = "This study"
cols = ["study", "Dot"]

combined = pd.concat([df.rename(columns={"scores": "Dot"})[cols], wang_db[cols], wang_db_tms[cols]])

make_simple_boxplot(combined, "study", "Dot", text_width=60).savefig("paper_plots/Fig9a_studies.png", bbox_inches='tight')


df_wang2020_db_with_classes = load_sdf_into_df("../data/reference/wang2020/3411.sdf")
wang_db_merged = df_wang2020_db_with_classes.set_index("inchikey").join(wang_db.set_index("Inchikey "))[["superclass","class","subclass","Dot"]].dropna()
order = [
    "Lipids and lipid-like molecules",
    "Organoheterocyclic compounds",
    "Organic oxygen compounds",
    "Organic acids and derivatives",
    "Organic nitrogen compounds",
    "Hydrocarbons"
]

make_simple_boxplot(wang_db_merged, x="superclass", y="Dot", order=order).savefig("paper_plots/Fig9b_wang2020.png", bbox_inches='tight')

df_wang2022_db_tms_with_classes = load_sdf_into_df("../data/reference/wang2022/3412.sdf")
wang_db_tms_merged = df_wang2022_db_tms_with_classes.set_index("inchikey").join(wang_db_tms.set_index("inchikey")["Dot"])[["superclass","class","subclass","Dot"]].dropna()
order = [
    "Benzenoids",
    "Lipids and lipid-like molecules",
    "Organoheterocyclic compounds",
    "Phenylpropanoids and polyketides",
    "Organic oxygen compounds",
    "Organic acids and derivatives",
    "Organic nitrogen compounds",
    "Organosulfur compounds",
    "Organometallic compounds"
]
make_simple_boxplot(wang_db_tms_merged, x="superclass", y="Dot", order = order).savefig("paper_plots/Fig9c_wang2022.png", bbox_inches='tight')