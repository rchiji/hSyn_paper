import os
import h5py
import joblib
from glob import glob
import gc
import numpy as np
import pandas as pd
import polars as pl
import scipy.sparse as sp
from pathlib import Path
from tqdm import tqdm
import random
from sklearn.decomposition import PCA
import umap
import matplotlib
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import seaborn as sns


# Aggregation_segmentation_outputs
## load data
anno = pl.read_csv("00_src/annotations_full.txt", separator="\t")
anno = anno.rename({anno.columns[0]: "Donor"})

donor_knee = (
    anno.filter(pl.col("Joint") == "Knee")
        .get_column("Donor")
        .to_list()
)

dfs = []
for f in glob("../GNN/data/SLICTile_Features_241125/*.txt"):
    donor = os.path.splitext(os.path.basename(f))[0]
    df = pl.read_csv(f, separator="\t")
    df = df.with_columns(pl.col("ID").cast(str))
    df = df.with_columns(pl.lit(donor).alias("Donor"))
    dfs.append(df)

## combine all data
all_data_combined = pl.concat(dfs)

all_data_combined_knee = all_data_combined.filter(pl.col("Donor").is_in(donor_knee))
print(all_data_combined_knee.shape)
# -> (11770418, 66)

joblib.dump(all_data_combined_knee, "03_python_outs/pkl/all_data_combined_knee.pkl")


# PCA - UMAP
## PCA
feature_columns = all_data_combined_knee.columns[1:65]
features = all_data_combined_knee.select(pl.col(feature_columns)).to_numpy()

pca = PCA(n_components=64, random_state=123)
pca_result = pca.fit_transform(features)

std_devs = np.sqrt(pca.explained_variance_)
prop_var = pca.explained_variance_ratio_ 
cum_prop = np.cumsum(pca.explained_variance_ratio_)
pca_summary = pd.DataFrame({
    "Standard Deviation": std_devs,
    "Proportion of Variance": prop_var,
    "Cumulative Proportion": cum_prop
}, index=[f"PC{i+1}" for i in range(len(std_devs))])

with pd.option_context("display.max_rows", None):
    print(pca_summary)

tile_ids = all_data_combined_knee[:, 0].to_numpy()
donors = all_data_combined_knee[:, 65].to_numpy()
pca_df = pd.DataFrame(pca_result, columns=[f"PC{i+1}" for i in range(pca_result.shape[1])])
pca_df["TileID"] = tile_ids
pca_df["Donor"] = donors

joblib.dump(pca_df, "03_python_outs/pkl/pca_df_20250713.pkl")
joblib.dump(pca_result, "03_python_outs/pkl/pca_model_20250713.pkl")

pca_result_select = pca_result[:, :46]


## UMAP
reducer = umap.UMAP(n_neighbors=15, min_dist=0.05, n_components=2, random_state=123)
umap_result = reducer.fit_transform(pca_result_select)

umap_df = pd.DataFrame(umap_result, columns=["UMAP1", "UMAP2"])
umap_df["TileID"] = tile_ids
umap_df["Donor"] = donors

joblib.dump(umap_df, "03_python_outs/pkl/umap_df_20250815.pkl")
joblib.dump(umap_result, "03_python_outs/pkl/umap_model_20250815.pkl")


# Visualization
df_all = umap_df.copy()
df_all["TileID"] = df_all["TileID"].astype(str)

label = []

for donor in df_all["Donor"].unique():
    with h5py.File("../GNN/GNN_InOut.h5", mode="r") as f:
        idx = f[donor]["index"][:].astype(str)
        maxFeature = f[donor]["tile_maxLabel"][:].astype(str)
    for i, lab in zip(idx, maxFeature):
        label.append({"Donor": donor, "ID": i, "Label": lab})
id_label_df = pd.DataFrame(label)

id_label_df["ID"] = id_label_df["ID"].astype(str)

df_all = df_all.merge(
    id_label_df[["ID", "Label"]],
    left_on="TileID", right_on="ID",
    how="left"
)

joblib.dump(df_all, "03_python_outs/pkl/combined_df_umap_label_leiden.pkl")


fm.fontManager.addfont("/usr/share/fonts/truetype/msttcorefonts/Arial.ttf")
fm.fontManager.addfont("/usr/share/fonts/truetype/msttcorefonts/Arial_Bold.ttf")
fm.fontManager.addfont("/usr/share/fonts/truetype/msttcorefonts/Arial_Italic.ttf")
fm.fontManager.addfont("/usr/share/fonts/truetype/msttcorefonts/Arial_Bold_Italic.ttf")

plt.rcParams.update({
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "font.family": "Arial",
    "font.sans-serif": ["Arial"],
})

rename_map = {
    "Fibro(dense,irregular)": "Fibrous tissue\n(dense, irregular)",
    "Fibro(dense,regular)":   "Fibrous tissue\n(dense, regular)",
    "Fibro(loose)":           "Fibrous tissue\n(loose)",
    "vessel":                 "Micro vessel",
    "vessel(large)":          "Large vessel",
    "lining":                 "Lining",
    "Immune cells":           "TLS",
    "plasma":                 "Plasma",
    "adipose":                "Adipose",
    "muscle":                 "Muscle",
    "Stroma":                 "Stroma",
    "RBC":                    "RBC"
}

label_order = [
    "Adipose",
    "Fibrous tissue\n(loose)",
    "Fibrous tissue\n(dense, irregular)",
    "Fibrous tissue\n(dense, regular)",
    "Stroma",
	"TLS",
    "Plasma",
    "Micro vessel",
    "Large vessel",	
	"RBC",
    "Lining",	
    "Muscle"
]

POINT_SIZE = 0.1
ALPHA = 0.8
LINEWIDTH_SPINE = 0.5

df_all_random = df_all.sample(n=500_000, random_state=42)
df_all_random["Label_renamed"] = df_all_random["Label"].replace(rename_map)

paired = sns.color_palette("Paired", n_colors=len(label_order))
palette_dict = dict(zip(label_order, paired))

df_plot = df_all_random.copy()
df_plot["Label_renamed"] = pd.Categorical(
    df_plot["Label_renamed"],
    categories=label_order,
    ordered=True
)

df_plot = df_plot.sort_values("Label_renamed")

fig, ax = plt.subplots(figsize=(2.5, 2.25))
fig.subplots_adjust(right=0.82)
sns.scatterplot(
    data=df_plot, x="UMAP1", y="UMAP2",
    hue="Label_renamed", hue_order=label_order, 
    palette=palette_dict,
    s=POINT_SIZE, linewidth=0, alpha=ALPHA, ax=ax, legend="full"
)

for collection in ax.collections:
    collection.set_rasterized(True)

ax.set_xlabel("UMAP1", fontsize=5)
ax.set_ylabel("UMAP2", fontsize=5)
ax.tick_params(axis="both", labelsize=4, width=0.4, length=2)
for spine in ax.spines.values():
    spine.set_visible(True)
    spine.set_linewidth(LINEWIDTH_SPINE)

ax.legend(
    bbox_to_anchor=(1.02, 1), loc="upper left",
    borderaxespad=0., title_fontsize=5, fontsize=5, title="Label",
    frameon=False, markerscale=5
)

plt.savefig("99_Fig/fig5/umap_by_label.png", dpi=300, bbox_inches="tight")
plt.savefig("99_Fig/fig5/umap_by_label.pdf", dpi=300, bbox_inches="tight")

