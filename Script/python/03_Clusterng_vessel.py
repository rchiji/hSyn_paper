import os
import gc
import json
import random
import joblib
import faiss
import h5py
import matplotlib
import matplotlib.font_manager as fm
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import polars as pl
import igraph as ig
import seaborn as sns
import leidenalg as la
import scipy.sparse as sp
from tqdm import tqdm
from glob import glob
from pathlib import Path
from collections import defaultdict

all_data_combined_knee = joblib.load("tmp/all_data_combined_knee.pkl")


with h5py.File("GNN_InOut.h5", "r") as f:
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

# joblib.dump(all_data_combined_knee, "pkl/all_data_combined_knee_with_label.pkl")


# Clustering
OUTDIR = "02_Clustering/leiden_pipeline_outs/vessel"
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


feature_cols = all_data_combined_knee_vessel.columns[1:65]
meta_cols = [c for c in all_data_combined_knee_vessel.columns if c not in feature_cols]
meta_df = all_data_combined_knee_vessel.select(meta_cols)

X = all_data_combined_knee_vessel.select(feature_cols).to_numpy()
if X.dtype != np.float32:
    X = X.astype(np.float32, copy=False)
X = np.ascontiguousarray(X)


if USE_COSINE and NORMALIZE_FOR_COSINE:
    norms = np.linalg.norm(X, axis=1, keepdims=True)
    X = X / (norms + 1e-12)
    print("Applied L2 normalization for cosine similarity.")

# np.save(os.path.join(OUTDIR, "normalized_features_vessel.npy"), X)


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

# np.save(os.path.join(OUTDIR, f"knn_indices_K{K}.npy"), I)
# np.save(os.path.join(OUTDIR, f"knn_sims_K{K}.npy"), Dsim)


## Check
# N, D = X.shape
# print(f"Feature matrix: {X.shape}, dtype={X.dtype}")
# print(I.shape, I.dtype, Dsim.shape, Dsim.dtype)
# 
# assert I.min() >= 0 and I.max() < N, "neighbor index out of range"
# self_hit = (I == np.arange(N, dtype=I.dtype)[:, None]).any()
# print("Has self neighbor?", bool(self_hit))  # False expected
# 
# assert np.isfinite(Dsim).all(), "NaN/Inf detected in similarities"
# print("sim min/max:", float(Dsim.min()), float(Dsim.max()))
# if MAKE_WEIGHTS_NONNEG or SHIFT_TO_01:
#     assert float(Dsim.min()) >= 0.0, "Weights must be non-negative for stability."
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

# np.save(os.path.join(OUTDIR, f"edges_src_K{K}.npy"), final_i)
# np.save(os.path.join(OUTDIR, f"edges_dst_K{K}.npy"), final_j)
# np.save(os.path.join(OUTDIR, f"edges_w_K{K}.npy"), final_w)


## Check
# print("Constructing igraph graph...")
# g = ig.Graph(n=N, directed=False)
# BATCH = 5_000_000
# for start in tqdm(range(0, len(final_i), BATCH)):
#     end = min(start + BATCH, len(final_i))
#     edges_batch = list(zip(final_i[start:end].tolist(), final_j[start:end].tolist()))
#     g.add_edges(edges_batch)
#     del edges_batch
# gc.collect()
# # weight
# assert g.ecount() == len(final_w), "Edge count mismatch for weight assignment."
# g.es["weight"] = final_w.tolist()
# print(g.summary())
# comp = g.components(mode="WEAK")
# gcc = comp.giant().vcount() / N
# print(f"GCC ratio: {gcc:.4%}")
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

# joblib.dump(partitions, os.path.join(OUTDIR, "leiden_res_comparison.pkl"))


# Visualization
RES = 0.05
OUT_CLUSTER_TAG = f"Vessel_Leiden_res{RES}" 
OUT_DIR_TXT = f"03_Community/SLICTile/SLICTile_cluster_res{RES}_vessel"
os.makedirs(OUT_DIR_TXT, exist_ok=True)

lab = pl.read_csv(f"02_Clustering/leiden_pipeline_outs/vessel/labels_with_meta_res{RES}.tsv", separator="\t")
assert "ID" in lab.columns and "Donor" in lab.columns, "TileID/Donor: not found"
cl_col = f"cluster_res_{RES}"
assert cl_col in lab.columns, f"{cl_col} : not found"
lab_df = lab.select(["ID", "Donor", cl_col]).to_pandas()
lab_df[cl_col] = lab_df[cl_col].astype(int)

## Reconstruct info_dict
info_dict = {}
with h5py.File("GNN_InOut.h5", mode="r") as f:
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

## Assign Leiden clusters by donor × TileID
for donor, d in tqdm(info_dict.items(), desc="Attach Leiden"):
    sub = lab_df[lab_df["Donor"] == donor].set_index("ID")
    cl = sub.reindex(d["index"])[cl_col]
    missing = cl.isna().sum()
    if missing > 0:
        print(f"[WARN] {donor}: {missing} tiles have no Leiden label (TileID: mismatch)")
    info_dict[donor][OUT_CLUSTER_TAG] = cl.fillna(-1).astype(int).values
valid_donors = set(lab_df["Donor"].unique())
info_dict = {d:obj for d, obj in info_dict.items() if d in valid_donors}

## Join cluster numbers, features, and input composition across all tiles
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

## Medoid extraction
medoids = {} 
for cl, ids in cluster_series.groupby(cluster_series).groups.items():
    if cl < 0:
        continue
    ids = list(ids)
    F = X_feat.loc[ids]
    centroid = F.mean(axis=0).values
    dist = np.linalg.norm(F.values - centroid[None, :], axis=1)
    argmin = np.argmin(dist)
    tile_id = F.index[argmin]
    donor = id_to_donor.at[tile_id]
    medoids[cl] = (donor, tile_id)

## Proportion of 'within-tile composition' for each cluster
row_sums = X_in.sum(axis=1).replace(0, np.nan)
X_in_norm = X_in.div(row_sums, axis=0).fillna(0.0)
df_ratio_merge_mean_tile = X_in_norm.groupby(cluster_series).mean() * 100
className = ["Immune cells", "plasma", "Fibro(loose)", "Fibro(dense,regular)",
             "Fibro(dense,irregular)", "lining", "vessel", "vessel(large)",
             "adipose", "Stroma", "muscle", "RBC"]
if df_ratio_merge_mean_tile.shape[1] == len(className):
    df_ratio_merge_mean_tile.columns = className

# Proportion of 'community (5-hop) composition' for each cluster
df_comm_ratio_list = []
for donor, d in info_dict.items():
    with h5py.File("GNN_InOut.h5", mode="r") as f:
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
    "Fibro(dense,irregular)":         "Fibrous tissue\n(dense, irregular)",
    "Fibro(dense,regular)":           "Fibrous tissue\n(dense, regular)",
    "Fibro(loose)":                   "Fibrous tissue\n(loose)",
    "vessel":                         "Micro vessel",
    "vessel(large)":                  "Large vessel",
    "lining":                         "Lining",
    "Immune cells":                   "TLS",
    "plasma":                         "Plasma",
    "adipose":                        "Adipose",
    "muscle":                         "Muscle"
}

def scale_rows_zscore(mat: pd.DataFrame, eps=1e-9):
    mean = mat.mean(axis=1)
    std  = mat.std(axis=1, ddof=1).replace(0, eps)
    return mat.sub(mean, axis=0).div(std, axis=0)

def heatmap_percent_modified_2(df, color, title_text, out_name, scale_row=True):
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

heatmap_percent_modified_2(
    df_ratio_merge_mean_tile,
    color="RdBu_r",
    title_text="Center tile",
    out_name=f"03_Community/Proportion/vessel/TileComp_{OUT_CLUSTER_TAG}.pdf",
	scale_row=True
)
heatmap_percent_modified_2(
    df_ratio_merge_mean_comm,
    color="RdBu_r",
    title_text="Tissue community (5-hop)",
    out_name=f"03_Community/Proportion/vessel/CommComp_{OUT_CLUSTER_TAG}.pdf",
	scale_row=True
)

## List of representative tiles (for preparation such as extraction in QuPath)
medoid_df = (pd.DataFrame.from_dict(medoids, orient="index", columns=["Donor","TileID"])
               .sort_index().rename_axis("Cluster").reset_index())
medoid_df.to_csv(f"leiden_pipeline_res/vessel/Medoids_{OUT_CLUSTER_TAG}.tsv", sep="\t", index=False)
print(medoid_df.head())

with h5py.File("GNN_InOut.h5", mode="a") as f:
    for donor, d in tqdm(info_dict.items(), desc="Write H5 Leiden"):
        ds_name = OUT_CLUSTER_TAG 
        if ds_name in f[donor]:
            del f[donor][ds_name]
        f[donor].create_dataset(ds_name, data=d[OUT_CLUSTER_TAG].astype(np.int32), compression="gzip")

for donor, d in info_dict.items():
    write_df = pd.DataFrame({"ID": d["index"], "Cluster": d[OUT_CLUSTER_TAG]})
    write_df.to_csv(os.path.join(OUT_DIR_TXT, f"{donor}.txt"), sep="\t", index=False)

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

out_dir = "03_Community/Map/Map_CellType_Leiden0.05_vessel"
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

order_cols = ["Immune cells","plasma","Fibro(loose)","Fibro(dense,regular)",
              "Fibro(dense,irregular)","lining","vessel","vessel(large)",
              "adipose","Stroma","muscle","RBC"]
M_tile = (df_ratio_merge_mean_tile
          .reindex(columns=order_cols)
          .sort_index())
M_tile.to_csv(f"leiden_pipeline_res/vessel/Matrix_TileComp_res{RES}_vessel.tsv", sep="\t", float_format="%.3f")
M_comm = (df_ratio_merge_mean_comm
          .reindex(columns=order_cols)
          .sort_index())
M_comm.to_csv(f"leiden_pipeline_res/vessel/Matrix_CommComp_res{RES}_vessel.tsv", sep="\t", float_format="%.3f")

