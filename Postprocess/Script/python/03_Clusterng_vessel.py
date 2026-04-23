import os
import h5py
import joblib
import gc
import json
import numpy as np
import pandas as pd
import polars as pl
import scipy.sparse as sp
from pathlib import Path
from tqdm import tqdm
import random
from sklearn.decomposition import PCA
import umap
import faiss
import igraph as ig
import leidenalg as la
import matplotlib
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import seaborn as sns


all_data_combined_knee = joblib.load("03_python_outs/pkl/all_data_combined_knee.pkl")


with h5py.File("../GNN/GNN_InOut.h5", "r") as f:
    dfs = []
    for case in f.keys(): 
        index = f[case]["index"][:]
        labels = f[case]["tile_maxLabel"][:].astype(str)
        df_case = pd.DataFrame({
            "index": index,
            "label": labels
        })
        dfs.append(df_case)   
    df_ID_label = pd.concat(dfs, axis=0, ignore_index=True)
keep = (all_data_combined_knee.select(pl.col("ID").unique()).to_series().to_list())
keep = set(x.encode("utf-8") for x in keep)
df_ID_label = df_ID_label[df_ID_label["index"].isin(keep)]
df_ID_label = (
    pl.from_pandas(df_ID_label)
    .with_columns(
        pl.col("index").cast(pl.Utf8, strict=False).alias("index_str"),
        pl.col("label").cast(pl.Utf8)
    )
)
all_data_combined_knee = all_data_combined_knee.join(
    df_ID_label.select(["index_str", "label"]),
    left_on="ID",
    right_on="index_str",
    how="left"
)

all_data_combined_knee_vessel = all_data_combined_knee.filter(
    pl.col("label").is_in(["vessel","vessel(large)"])
)

joblib.dump(all_data_combined_knee, "03_python_outs/pkl/all_data_combined_knee_with_label.pkl")



# Dimensionality_reduction
feature_columns = all_data_combined_knee_vessel.columns[1:65]
features = all_data_combined_knee_vessel.select(pl.col(feature_columns)).to_numpy()

pca = PCA(n_components=64, random_state=123)
pca_result = pca.fit_transform(features)

plt.plot(np.arange(1, 65), pca.explained_variance_, marker='o')
plt.title("PCA - Explained Variance per Component")
plt.xlabel("Principal Component")
plt.ylabel("Explained Variance")
plt.grid(True)
plt.show()

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

tile_ids = all_data_combined_knee_vessel[:, 0].to_numpy()
donors = all_data_combined_knee_vessel[:, 65].to_numpy()
pca_df = pd.DataFrame(pca_result, columns=[f"PC{i+1}" for i in range(pca_result.shape[1])])
pca_df["TileID"] = tile_ids
pca_df["Donor"] = donors

joblib.dump(pca_df, "03_python_outs/pkl/pca_df_vessel_20260318.pkl")
joblib.dump(pca_result, "03_python_outs/pkl/pca_model_vessel_20260318.pkl")

pca_result_select = pca_result[:, :13]

reducer = umap.UMAP(n_neighbors=15, min_dist=0.05, n_components=2, random_state=123)
umap_result = reducer.fit_transform(pca_result_select)
umap_df = pd.DataFrame(umap_result, columns=["UMAP1", "UMAP2"])
umap_df["TileID"] = tile_ids
umap_df["Donor"] = donors
joblib.dump(umap_df, "03_python_outs/pkl/umap_df_vessel_20260318.pkl")
joblib.dump(umap_result, "03_python_outs/pkl/umap_model_vessel_20260318.pkl")



# Clustering
OUTDIR = "03_python_outs/Clustering/leiden_pipeline_outs/vessel"
os.makedirs(OUTDIR, exist_ok=True)
## neighbor search
K = 20
USE_COSINE = True
NORMALIZE_FOR_COSINE = True 
HNSW_M = 32
HNSW_EF_CONSTRUCT = 200
HNSW_EF_SEARCH = 128
## clustering
RESOLUTIONS = [0.01, 0.05, 0.075, 0.1, 0.2]
SEED = 123
random.seed(SEED)
np.random.seed(SEED)
MAKE_WEIGHTS_NONNEG = True
SHIFT_TO_01 = False
SYMMETRIZE_MODE = "union"


feature_columns = all_data_combined_knee_vessel.columns[1:65]
meta_cols = [c for c in all_data_combined_knee_vessel.columns if c not in feature_columns]
meta_df = all_data_combined_knee_vessel.select(meta_cols)

X = all_data_combined_knee_vessel.select(feature_columns).to_numpy()
if X.dtype != np.float32:
    X = X.astype(np.float32, copy=False)

X = np.ascontiguousarray(X)


if USE_COSINE and NORMALIZE_FOR_COSINE:
    norms = np.linalg.norm(X, axis=1, keepdims=True)
    X = X / (norms + 1e-12)

np.save(os.path.join(OUTDIR, "normalized_features_vessel.npy"), X)


metric = faiss.METRIC_INNER_PRODUCT if USE_COSINE else faiss.METRIC_L2
index = faiss.IndexHNSWFlat(X.shape[1], HNSW_M, metric)
index.hnsw.efConstruction = HNSW_EF_CONSTRUCT
if hasattr(index.hnsw, "seed"):
    index.hnsw.seed = SEED

index.add(X)
index.hnsw.efSearch = HNSW_EF_SEARCH

distances, indices = index.search(X, K + 1)
I = indices[:, 1:]
Dsim = distances[:, 1:]
if USE_COSINE:
    if MAKE_WEIGHTS_NONNEG and not SHIFT_TO_01:
        np.maximum(Dsim, 0.0, out=Dsim)
    elif SHIFT_TO_01:
        Dsim = (Dsim + 1.0) * 0.5  
else:
    Dsim = 1.0 / (1.0 + Dsim) 

np.save(os.path.join(OUTDIR, f"knn_indices_K{K}.npy"), I)
np.save(os.path.join(OUTDIR, f"knn_sims_K{K}.npy"), Dsim)


## Check
N, D = X.shape
print(f"Feature matrix: {X.shape}, dtype={X.dtype}")
print(I.shape, I.dtype, Dsim.shape, Dsim.dtype)

assert I.min() >= 0 and I.max() < N, "neighbor index out of range"
self_hit = (I == np.arange(N, dtype=I.dtype)[:, None]).any()
print("Has self neighbor?", bool(self_hit))  # False expected

assert np.isfinite(Dsim).all(), "NaN/Inf detected in similarities"
print("sim min/max:", float(Dsim.min()), float(Dsim.max()))
if MAKE_WEIGHTS_NONNEG or SHIFT_TO_01:
    assert float(Dsim.min()) >= 0.0, "Weights must be non-negative for stability."
# -> Feature matrix: (201631, 64), dtype=float32
# -> (201631, 20) int64 (201631, 20) float32
# -> Has self neighbor? False
# -> sim min/max: 0.07825250178575516 0.9991085529327393


row_ids = np.repeat(np.arange(N, dtype=np.int64), K)
col_ids = I.reshape(-1).astype(np.int64)
w_vals  = Dsim.reshape(-1).astype(np.float32)

mask = row_ids != col_ids
row_ids = row_ids[mask]; col_ids = col_ids[mask]; w_vals = w_vals[mask]

W = sp.coo_matrix((w_vals, (row_ids, col_ids)), shape=(N, N), dtype=np.float32)

if SYMMETRIZE_MODE == "mutual":
    W_sym = (W.minimum(W.T)).tocsr()  # Mutual neighbors only
else:
    W_sym = (W.maximum(W.T)).tocsr() 
W_sym.setdiag(0.0)
W_sym.eliminate_zeros()

W_coo = W_sym.tocoo()
upper_mask = W_coo.row < W_coo.col
final_i = W_coo.row[upper_mask].astype(np.int32)
final_j = W_coo.col[upper_mask].astype(np.int32)
final_w = W_coo.data[upper_mask].astype(np.float32)
print(f"Edges (unique, undirected): {final_i.size:,}")
assert np.all(final_w >= 0.0), "Negative weights detected after symmetrization."

np.save(os.path.join(OUTDIR, f"edges_src_K{K}.npy"), final_i)
np.save(os.path.join(OUTDIR, f"edges_dst_K{K}.npy"), final_j)
np.save(os.path.join(OUTDIR, f"edges_w_K{K}.npy"), final_w)


## Check
g = ig.Graph(n=N, directed=False)
BATCH = 5_000_000
for start in tqdm(range(0, len(final_i), BATCH)):
    end = min(start + BATCH, len(final_i))
    edges_batch = list(zip(final_i[start:end].tolist(), final_j[start:end].tolist()))
    g.add_edges(edges_batch)
    del edges_batch
gc.collect()
### weight
assert g.ecount() == len(final_w), "Edge count mismatch for weight assignment."
g.es["weight"] = final_w.tolist()
print(g.summary())
comp = g.components(mode="WEAK")
gcc = comp.giant().vcount() / N
print(f"GCC ratio: {gcc:.4%}")
# -> Constructing igraph graph...
# -> 100%|██████████| 1/1 [00:00<00:00,  1.51it/s]
# -> IGRAPH U-W- 201631 3119870 -- 
# -> + attr: weight (e)
# -> GCC ratio: 100.0000%


partitions = {}
for res in RESOLUTIONS:
    print(f"Running Leiden (RBConfigurationVertexPartition, resolution={res}) ...")
    partition_type = la.RBConfigurationVertexPartition  
    part = la.find_partition(
        g,
        partition_type,
        weights=g.es["weight"],
        resolution_parameter=res,
        seed=SEED
    )
    labels = np.array(part.membership, dtype=np.int32)
    np.save(os.path.join(OUTDIR, f"leiden_labels_res{res}.npy"), labels)
    # summary
    n_clusters = int(np.unique(labels).size)
    sizes = np.bincount(labels)
    print(f"resolution={res}: clusters={n_clusters}, min={sizes.min()}, median={np.median(sizes)}, max={sizes.max()}")
    partitions[res] = {
        "n_clusters": n_clusters,
        "min_size": int(sizes.min()),
        "median_size": float(np.median(sizes)),
        "max_size": int(sizes.max()),
    }

with open(os.path.join(OUTDIR, "leiden_partitions_summary.json"), "w") as f:
    json.dump(partitions, f, indent=2)

for res in RESOLUTIONS:
    labels = np.load(os.path.join(OUTDIR, f"leiden_labels_res{res}.npy"))
    out_pl = meta_df.with_columns(pl.Series(name=f"cluster_res_{res}", values=labels))
    out_pl.write_csv(os.path.join(OUTDIR, f"labels_with_meta_res{res}.tsv"), separator="\t")

joblib.dump(partitions, os.path.join(OUTDIR, "leiden_res_comparison.pkl"))


# Visualization
## UMAP
LABELDIR = Path("03_python_outs/Clustering/leiden_pipeline_outs/vessel")
resolutions = [0.05, 0.1]

df_all = umap_df.copy()
df_all["TileID"] = df_all["TileID"].astype(str)

for res in resolutions:
    lab = pl.read_csv(LABELDIR / f"labels_with_meta_res{res}.tsv", separator="\t")
    lab_df = lab.to_pandas()
    lab_df["ID"] = lab_df["ID"].astype(str)
    colname = f"cluster_res_{res}"
    df_all = df_all.merge(lab_df[["ID", colname]], left_on="TileID", right_on="ID", how="left")
    df_all.drop(columns=["ID"], inplace=True)

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

joblib.dump(df_all, "03_python_outs/pkl/combined_df_umap_label_leiden_vessel.pkl")


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

df_all["Label_renamed"] = df_all["Label"].replace(rename_map)

POINT_SIZE = 0.1
ALPHA = 0.8
LINEWIDTH_SPINE = 0.5

col = f"cluster_res_0.05"
turbo3_hex = ["#30123B", "#FABA39", "#7A0403"]
palette = {0: turbo3_hex[0], 1: turbo3_hex[1], 2: turbo3_hex[2]}

fig, ax = plt.subplots(figsize=(2.5, 2.25))
fig.subplots_adjust(right=0.82)
sns.scatterplot(
    data=df_all, x="UMAP1", y="UMAP2",
    hue=col, palette=palette,
    s=POINT_SIZE, linewidth=0, alpha=ALPHA, ax=ax, legend="full"
)

for collection in ax.collections:
    collection.set_rasterized(True)

for side in ["top", "right", "bottom", "left"]:
    ax.spines[side].set_visible(True)
    ax.spines[side].set_linewidth(LINEWIDTH_SPINE)

ax.set_xlabel("UMAP1", fontsize=5)
ax.set_ylabel("UMAP2", fontsize=5)
ax.tick_params(axis="both", labelsize=4, width=0.4, length=2)

ax.legend(
    bbox_to_anchor=(1.02, 1), loc="upper left",
    borderaxespad=0., title_fontsize=5, fontsize=5, title="Cluster",
    frameon=False, markerscale=5
)

out_png = f"99_Fig/fig5/umap_res0.05_leiden_random_vessel.png"
out_pdf = f"99_Fig/fig5/umap_res0.05_leiden_random_vessel.pdf"

fig.savefig(out_png, dpi=300, bbox_inches="tight")
fig.savefig(out_pdf, dpi=300, bbox_inches="tight")

plt.close(fig)

paired = sns.color_palette("Paired", n_colors=len(label_order))
palette_dict = dict(zip(label_order, paired))

df_plot = df_all.copy()
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

plt.savefig("99_Fig/fig5/umap_by_label_vessel.png", dpi=300, bbox_inches="tight")
plt.savefig("99_Fig/fig5/umap_by_label_vessel.pdf", dpi=300, bbox_inches="tight")
plt.close(fig)


## Community
RES = 0.05
OUT_CLUSTER_TAG = f"Vessel_Leiden_res{RES}" 
OUT_DIR_TXT = f"03_python_outs/Community/SLICTile/SLICTile_cluster_res{RES}_vessel"
os.makedirs(OUT_DIR_TXT, exist_ok=True)

lab = pl.read_csv(f"03_python_outs/Clustering/leiden_pipeline_outs/vessel/labels_with_meta_res{RES}.tsv", separator="\t")
assert "ID" in lab.columns and "Donor" in lab.columns, "TileID/Donor: not found"
cl_col = f"cluster_res_{RES}"
assert cl_col in lab.columns, f"{cl_col} : not found"
lab_df = lab.select(["ID", "Donor", cl_col]).to_pandas()
lab_df[cl_col] = lab_df[cl_col].astype(int)

### Reconstruct info_dict
info_dict = {}
with h5py.File("../GNN/GNN_InOut.h5", mode="r") as f:
    donors = list(f.keys())
    for donor in tqdm(donors, desc="Load H5"):
        idx = f[donor]["index"][:].astype(str)
        coords = f[donor]["coord"][:]
        input_df = pd.DataFrame(f[donor]["input_feat"][:], index=idx)
        features = pd.DataFrame(f[donor]["pred_features"][:], index=idx)
        maxID = f[donor]["tile_maxID"][:]
        maxFeature = f[donor]["tile_maxLabel"][:].astype(str)
        mask = (maxFeature == "vessel") | (maxFeature == "vessel(large)")
        idx = idx[mask]
        coords = coords[mask]
        input_df = input_df.loc[idx]
        features = features.loc[idx]
        maxID = maxID[mask]
        maxFeature = maxFeature[mask]
        info_dict[donor] = {
            "index": idx,
            "coord": pd.DataFrame(coords, index=idx, columns=["x","y"]),
            "input_df": input_df,
            "features": features,
            "id_max": maxID,
        }

### Assign Leiden clusters by donor × TileID
for donor, d in tqdm(info_dict.items(), desc="Attach Leiden"):
    sub = lab_df[lab_df["Donor"] == donor].set_index("ID")
    cl = sub.reindex(d["index"])[cl_col]
    missing = cl.isna().sum()
    if missing > 0:
        print(f"[WARN] {donor}: {missing} tiles have no Leiden label (TileID: mismatch)")
    info_dict[donor][OUT_CLUSTER_TAG] = cl.fillna(-1).astype(int).values
valid_donors = set(lab_df["Donor"].unique())
info_dict = {d:obj for d, obj in info_dict.items() if d in valid_donors}

### Join cluster numbers, features, and input composition across all tiles
all_features = []
all_input = []
all_idx = []
all_clusters = []
for donor, d in info_dict.items():
    idx = d["index"]
    all_idx.append(pd.DataFrame({"Donor": donor, "ID": idx}))
    all_features.append(pd.DataFrame(d["features"].values, index=idx))
    all_input.append(pd.DataFrame(d["input_df"].values, index=idx))
    all_clusters.append(pd.Series(d[OUT_CLUSTER_TAG], index=idx, name="Cluster"))
df_idx = pd.concat(all_idx, ignore_index=True) 
X_feat = pd.concat(all_features, axis=0)
X_in = pd.concat(all_input, axis=0)
cluster_series = pd.concat(all_clusters, axis=0)
assert len(df_idx) == len(X_feat) == len(X_in) == len(cluster_series)
common_idx = cluster_series.index.intersection(X_feat.index).intersection(X_in.index)
cluster_series = cluster_series.loc[common_idx]
X_feat = X_feat.loc[common_idx]
X_in   = X_in.loc[common_idx]
id_to_donor = df_idx.set_index("ID")["Donor"]
assert id_to_donor.index.is_unique, "ID: not unique"
assert set(cluster_series.index).issubset(set(id_to_donor.index)), "ID→Donor: deficient"

### Proportion of 'within-tile composition' for each cluster
row_sums = X_in.sum(axis=1).replace(0, np.nan)
X_in_norm = X_in.div(row_sums, axis=0).fillna(0.0)
df_ratio_merge_mean_tile = X_in_norm.groupby(cluster_series).mean() * 100
className = ["Immune cells", "plasma", "Fibro(loose)", "Fibro(dense,regular)",
             "Fibro(dense,irregular)", "lining", "vessel", "vessel(large)",
             "adipose", "Stroma", "muscle", "RBC"]
if df_ratio_merge_mean_tile.shape[1] == len(className):
    df_ratio_merge_mean_tile.columns = className

### Proportion of 'community (5-hop) composition' for each cluster
df_comm_ratio_list = []
for donor, d in info_dict.items():
    with h5py.File("../GNN/GNN_InOut.h5", mode="r") as f:
        ratio = f[donor]["community_label_ratio"][:]
        idx = f[donor]["index"][:].astype(str)
        maxFeature = f[donor]["tile_maxLabel"][:].astype(str)
        mask = (maxFeature == "vessel") | (maxFeature == "vessel(large)")
        idx = idx[mask]
        ratio = ratio[mask]
        maxFeature = maxFeature[mask]
    df = pd.DataFrame(ratio, index=idx)
    rs = df.sum(axis=1).replace(0, np.nan)
    df_norm = df.div(rs, axis=0).fillna(0.0)
    df_norm["Cluster"] = d[OUT_CLUSTER_TAG]
    df_comm_ratio_list.append(df_norm)

df_comm_ratio = pd.concat(df_comm_ratio_list, axis=0)
df_ratio_merge_mean_comm = df_comm_ratio.groupby("Cluster").mean().drop(columns=[], errors='ignore') * 100
if df_ratio_merge_mean_comm.shape[1] == len(className):
    df_ratio_merge_mean_comm.columns = className


def scale_rows_zscore(mat: pd.DataFrame, eps=1e-9):
    mean = mat.mean(axis=1)
    std  = mat.std(axis=1, ddof=1).replace(0, eps)
    return mat.sub(mean, axis=0).div(std, axis=0)

def rename_tissues_columns(df):
    df2 = df.copy()
    df2.columns = [rename_map.get(c, c) for c in df2.columns]
    return df2

def heatmap_percent_modified(df, color, title_text, out_name, scale_row=True):
    df_rename = rename_tissues_columns(df)
    mat = df_rename.T
    if scale_row:
        mat_plot = scale_rows_zscore(mat)
        vmin, vmax = -2, 2
        fmt = ".2f"
        cbar_label = "Row z-score"
    else:
        mat_plot = mat
        vmin, vmax = 0, 100
        fmt = ".1f"
        cbar_label = ""
    plt.figure(figsize=(6, 4))
    p = sns.heatmap(
        mat_plot,
        cmap=color,
        vmin=vmin, vmax=vmax,
        annot=True,
        fmt=fmt,
        annot_kws={"size": 4},
        linewidths=0.4, linecolor="#cfcfcf",
        cbar_kws={"shrink": 0.7},
        square=True
    )
    p.set_title(title_text, fontsize=6, pad=4)
    p.set_xlabel("Cluster", fontsize=5)
    p.set_ylabel("Tissue",  fontsize=5)
    p.tick_params(axis='x', labelsize=5, rotation=0, pad=1)
    p.tick_params(axis='y', labelsize=5, pad=1)
    if p.collections and p.collections[0].colorbar is not None:
        p.collections[0].colorbar.set_label(cbar_label, fontsize=5)
        p.collections[0].colorbar.ax.tick_params(labelsize=5)
    plt.tight_layout(pad=0.5)
    os.makedirs(os.path.dirname(out_name), exist_ok=True)
    plt.savefig(out_name)
    plt.show()

heatmap_percent_modified(
    df_ratio_merge_mean_tile,
    color="RdBu_r",
    title_text="Center tile",
    out_name=f"99_Fig/fig5/TileComp_{OUT_CLUSTER_TAG}.pdf",
	scale_row=True
)
heatmap_percent_modified(
    df_ratio_merge_mean_comm,
    color="RdBu_r",
    title_text="Tissue community (5-hop)",
    out_name=f"99_Fig/fig5/CommComp_{OUT_CLUSTER_TAG}.pdf",
	scale_row=True
)

for donor, d in info_dict.items():
    write_df = pd.DataFrame({"ID": d["index"], "Cluster": d[OUT_CLUSTER_TAG]})
    write_df.to_csv(os.path.join(OUT_DIR_TXT, f"{donor}.txt"), sep="\t", index=False)


## Spatial images
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

paired = sns.color_palette("Paired", n_colors=len(label_order))
palette_dict = dict(zip(label_order, paired))

turbo = matplotlib.cm.get_cmap("turbo")
turbo3_hex = ["#30123B", "#FABA39", "#7A0403"]

def plot_cluster(sample_info_dict, cluster_array, num_cluster,
                  pt_size=1, figsize=(14,5),
                  xmin=None, xmax=None, ymin=None, ymax=None,
                  className=None, color_map=None, cluster_cmap="tab20",
                  title_suffix="", flip=False, show=True):
    df = sample_info_dict["coord"]
    _label = sample_info_dict["id_max"]
    _cluster = np.asarray(cluster_array)
    xmin = df["x"].min() if xmin is None else xmin
    xmax = df["x"].max() if xmax is None else xmax
    ymin = df["y"].min() if ymin is None else ymin
    ymax = df["y"].max() if ymax is None else ymax
    _df = df.reset_index()
    ind = _df.query("(@xmin <= x <= @xmax) and (@ymin <= y <= @ymax)").index
    _df = _df.loc[ind]
    _cluster = _cluster[ind]
    _label = _label[ind]
    if flip:
        cy = (ymin + ymax) / 2
        _df["y"] = 2*cy - _df["y"]
    if className is None:
        className = list(range(1, int(np.max(_label))+1))
    display_names = [rename_map.get(name, name) for name in className]
    fig = plt.figure(figsize=figsize)
    ax0 = plt.subplot(121)
    for i, name in enumerate(className, start=1):
        idx = (_label == i)
        if not np.any(idx):
            continue
        label_name = display_names[i-1]
        if color_map is not None and label_name in color_map:
            point_color = color_map[label_name]
        else:
            point_color = "gray"
        ax0.scatter(
            _df.loc[idx, "x"], _df.loc[idx, "y"],
            s=pt_size, color=point_color,
            label=label_name
        )
    for col in ax0.collections:
        col.set_rasterized(True)
    ax0.set_title("Tissue Type" + title_suffix)
    handles, labels = ax0.get_legend_handles_labels()
    label_to_handle = dict(zip(labels, handles))
    ordered_labels = [lab for lab in label_order if lab in label_to_handle]
    ordered_handles = [label_to_handle[lab] for lab in ordered_labels]
    ax0.legend(
        ordered_handles, ordered_labels,
        title="",
        bbox_to_anchor=(1.01, 1), loc="upper left",
        borderaxespad=0., markerscale=4,
        fontsize=6, title_fontsize=6,
        frameon=False
    )
    ax1 = plt.subplot(122)
    cluster_color_map = {0: turbo3_hex[0], 1: turbo3_hex[1], 2: turbo3_hex[2]}

    for cid in np.sort(np.unique(_cluster)):
        idx = (_cluster == cid)
        ax1.scatter(
            _df.loc[idx, "x"], _df.loc[idx, "y"],
            c=cluster_color_map.get(int(cid), "#999999"),
            s=pt_size,
            label=f"Cluster {cid}"
        )
    ax1.set_title("Leiden clustering" + title_suffix)
    ax1.legend(
        title="",
        bbox_to_anchor=(1.01, 1), loc="upper left",
        borderaxespad=0., markerscale=4,
        fontsize=6, title_fontsize=6,
        frameon=False
    )
    for ax in (ax0, ax1):
        xlim = ax.get_xlim()
        ylim = ax.get_ylim()
        ax.set_aspect('equal', adjustable='box')
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)
    if show:
        plt.show()
    return fig

out_dir = "99_Fig/fig5/Map_CellType_Leiden0.05_vessel"
os.makedirs(out_dir, exist_ok=True)
class_names = [
    "Immune cells","plasma","Fibro(loose)","Fibro(dense,regular)",
    "Fibro(dense,irregular)","lining","vessel","vessel(large)",
    "adipose","Stroma","muscle","RBC"
]

for donor in sorted(valid_donors):
    if donor not in info_dict:
        continue
    cluster_arr = info_dict[donor][OUT_CLUSTER_TAG]
    n_clusters = int(np.unique(cluster_arr[cluster_arr >= 0]).size)
    fig = plot_cluster(
        sample_info_dict=info_dict[donor],
        cluster_array=cluster_arr,
        num_cluster=n_clusters,
        color_map=palette_dict, 
        className=class_names,
        cluster_cmap="tab20",
        title_suffix=f"({OUT_CLUSTER_TAG})",
        flip=True,
        show=False
    )
    out_path = os.path.join(out_dir, f"{donor}_{OUT_CLUSTER_TAG}.pdf")
    fig.savefig(out_path, bbox_inches="tight", dpi=300)
    plt.close(fig)

