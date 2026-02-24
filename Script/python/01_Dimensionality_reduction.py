import os
import gc
import json
import random
import joblib
import faiss
import h5py
import umap
import matplotlib
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import polars as pl
import igraph as ig
import seaborn as sns
import scipy.sparse as sp
from glob import glob
from tqdm import tqdm
from collections import defaultdict
from sklearn.decomposition import PCA


# Aggregation_segmentation_outputs
## load data
anno = pl.read_csv("data/annotations_full.txt", separator="\t")
anno = anno.rename({anno.columns[0]: "Donor"})

donor_knee = (
    anno.filter(pl.col("Joint") == "Knee")
        .get_column("Donor")
        .to_list()
)

dfs = []
for f in glob("data/SLICTile_Features_241125/*.txt"):
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

# joblib.dump(all_data_combined_knee, "pkl/all_data_combined_knee.pkl")


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

# joblib.dump(pca_df, "pkl/pca_df_20250713.pkl")
# joblib.dump(pca_result, "pkl/pca_model_20250713.pkl")

pca_result_select = pca_result[:, :46]


## UMAP
reducer = umap.UMAP(n_neighbors=15, min_dist=0.05, n_components=2, random_state=123)
umap_result = reducer.fit_transform(pca_result_select)

umap_df = pd.DataFrame(umap_result, columns=["UMAP1", "UMAP2"])
umap_df["TileID"] = tile_ids
umap_df["Donor"] = donors

# joblib.dump(umap_df, "pkl/umap_df_20250815.pkl")
# joblib.dump(umap_result, "pkl/umap_model_20250815.pkl")

