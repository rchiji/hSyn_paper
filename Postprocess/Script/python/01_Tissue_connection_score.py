import os
import pandas as pd
import numpy as np
import seaborn as sns
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from sklearn.neighbors import radius_neighbors_graph
from glob import glob
from tqdm import tqdm
from collections import Counter
from scipy.stats import gaussian_kde
from numpy import dot
from numpy.linalg import norm


filepaths = glob("../GNN/tile_label_coord_241115/*txt")

df_dict = {}

for filepath in filepaths:
    samplename = os.path.basename(filepath)[:-7]
    df_dict[samplename] = pd.read_table(filepath, index_col=0)


# Output path for edge information
edge_df_dict = {}

# Distance threshold to define a connection
radius = 80

for name, _df in df_dict.items():
    # Extract XY coordinates
    points = _df[["x","y"]]

    adj_matrix = radius_neighbors_graph(X=points, radius=radius)
    edge_index = np.argwhere(adj_matrix > 0)

    # Format into a DataFrame
    edge_df = pd.DataFrame(edge_index, columns=["edge_from","edge_to"])
    edge_df["id_from"] = _df.index[edge_df["edge_from"]]
    edge_df["id_to"] = _df.index[edge_df["edge_to"]]
    edge_df["label_from"] = _df.iloc[edge_df["edge_from"], _df.columns.get_loc("label") ].values
    edge_df["label_to"] = _df.iloc[edge_df["edge_to"], _df.columns.get_loc("label")].values

    edge_df_dict[name] = edge_df


with PdfPages('../GNN/radius_graph_241115.pdf') as pdf:

    for donor, _edge_df in tqdm(edge_df_dict.items()):
        _edge_index = _edge_df[["edge_from","edge_to"]]
        _df = df_dict[donor]
        _points = _df[["x","y"]].copy()
        _points["y"] = -_points["y"]
        plt.figure(figsize=(50,40))
        graph = nx.from_edgelist(_edge_index.values)
        nx.draw(graph, pos=_points.to_numpy(), node_size=1)
        plt.title(donor)
        # plt.tight_layout()
        pdf.savefig(bbox_inches="tight")
        plt.close()
plt.show()


def get_neighbors(
    center_node, edge_index=None, adj_marix=None
):
    if edge_index is not None:
        if edge_index.shape[0] < edge_index.shape[1]:
            edge_index = edge_index.T

        neighbor_nodes = edge_index[edge_index[:,0] == center_node][:,1]

    if adj_marix is not None:
        neighbor_nodes = np.where(adj_matrix[center_node,:])[0]

    return neighbor_nodes


dirname = "tile_info_241115"
os.makedirs(dirname, exist_ok=True)


cell_type_levels = ["Immune cells", "plasma", "Fibro(loose)", "Fibro(dense,regular)",
"Fibro(dense,irregular)", "Fibro(Perivascular)", "lining", "vessel", "vessel(large)",
"adipose", "nerve", "Stroma", "muscle", "RBC"]

for donor, _df in tqdm(df_dict.items()) :

    num_nodes = _df.shape[0]
    edge_df = edge_df_dict[donor]
    edge_index = edge_df[["edge_from","edge_to"]].to_numpy()
    degrees = []
    purities = []
    entropies = []

    for i in range(num_nodes):
        center_celltype = _df.iloc[i]["label"]
        neighbors = get_neighbors(center_node=i, edge_index=edge_index)
        degrees.append(len(neighbors))

        if len(neighbors) > 0:
            count_dict = { celltype:0 for celltype in cell_type_levels }
            count_dict[center_celltype] += 1
            count_dict.update( Counter(_df.iloc[neighbor]["label"] for neighbor in neighbors) )
        
            mode_label = max(count_dict, key=count_dict.get)
            purity = count_dict[mode_label] / sum(count_dict.values())

            total_count = sum(count_dict.values())
            probabilities = [ count / total_count for count in count_dict.values()]

            entropy = -sum([p * np.log(p) for p in probabilities if p > 0])

        else:
            purity = np.nan
            entropy = np.nan

        purities.append(purity)
        entropies.append(entropy)

    _df["Degree"] = degrees
    _df["Purity"] = purities
    _df["Entropy"] = entropies

    _df.to_csv(f"{dirname}/{donor}.txt",sep="\t",index_label="ID")


def mean_filter(
    metrics,
    num_nodes=None,
    edge_index=None,
    adj_matrix=None
):
    if isinstance(metrics,list):
        metrics = np.array(metrics)

    num_nodes = metrics.shape[0]

    mean_values = []
    for i in range(num_nodes):
        neighbors = get_neighbors(i,edge_index,adj_matrix)
        center_value = metrics[i]

        if len(neighbors) > 0:
            neighbor_values = metrics[neighbors]
            mean = np.mean(np.append(neighbor_values,center_value))
        else:
            mean = center_value

        mean_values.append(mean)

    return mean_values


dirname = "../GNN/tile_info_241115"
os.makedirs(dirname, exist_ok=True)


cell_type_levels = ["Immune cells", "plasma", "Fibro(loose)", "Fibro(dense,regular)",
"Fibro(dense,irregular)", "Fibro(Perivascular)", "lining", "vessel", "vessel(large)",
"adipose", "nerve", "Stroma", "muscle", "RBC"]

for donor, _df in tqdm(df_dict.items()) :

    num_nodes = _df.shape[0]
    edge_df = edge_df_dict[donor]
    edge_index = edge_df[["edge_from","edge_to"]].to_numpy()

    _df["Purity_mean"] = mean_filter(_df["Purity"].values, edge_index=edge_index)
    _df["Entropy_mean"] = mean_filter(_df["Entropy"].values, edge_index=edge_index)

    _df.to_csv(f"{dirname}/{donor}.txt",sep="\t",index_label="ID")


plt.hist( df_dict["D001"]["Entropy_mean"])


for donor, _df in df_dict.items():
    data = _df["Entropy_mean"]
    data = data.dropna()
    density =gaussian_kde(data)
    x = np.linspace(data.min(), data.max(), 1000)
    plt.plot(x, density(x), label='Density curve')

for donor, _df in df_dict.items():
    data = _df["Purity_mean"]
    data = data.dropna()
    density =gaussian_kde(data)
    x = np.linspace(data.min(), data.max(), 1000)
    plt.plot(x, density(x), label='Density curve')


x_grid = np.linspace(0,2,100)

densities = []
for donor, _df in df_dict.items():
    data = _df["Entropy_mean"]
    data = data.dropna()
    density = gaussian_kde(data)
    densities.append(density(x_grid))


n = len(densities)

cosine_matrix = np.zeros((n, n))

for i in range(n):
    for j in range(n):
        cosine_similarity = dot(densities[i], densities[j]) / (norm(densities[i]) * norm(densities[j]))
        cosine_matrix[i, j] = cosine_similarity 


cosine_matrix = pd.DataFrame(cosine_matrix, index=df_dict.keys(), columns=df_dict.keys())

sns.set(font_scale=0.6)
p = sns.clustermap(cosine_matrix, figsize=(10,10),xticklabels=1,yticklabels=1, dendrogram_ratio=(0.1,0.1),)
p.savefig("99_Fig/sup_fig6/001_Entropy_density_cosinesimilality_241115.pdf", format="pdf")


cosine_matrix.index[p.dendrogram_col.reordered_ind]
# Index(['D031', 'D108', 'D028', 'D059', 'D006', 'D093_L', 'D027', 'D061',
#        'D024', 'D075',
#        ...
#        'D002', 'D036', 'D106', 'D034', 'D128', 'D064', 'D100', 'D033', 'D045',
#        'D104'],
#       dtype='object', length=118)



output_pdf = "99_Fig/sup_fig6/002_density_cluster1.pdf"

donors = cosine_matrix.index[p.dendrogram_col.reordered_ind][0:16]

with PdfPages(output_pdf) as pdf:
    plt.figure(figsize=(10, 6))

    for donor, _df in df_dict.items():
        if donor in donors:
            data = _df["Entropy_mean"].dropna()
            density = gaussian_kde(data)
            x = np.linspace(data.min(), data.max(), 100)
            plt.plot(x, density(x), label=f'{donor}')

    plt.xlabel("X-axis Label")
    plt.ylabel("Density")
    plt.title("Overlayed Density Plots")

    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small', ncol=1)
    plt.grid(True)

    pdf.savefig(bbox_inches='tight')
    plt.close()


output_pdf = "99_Fig/sup_fig6/002_density_cluster2.pdf"

donors = cosine_matrix.index[p.dendrogram_col.reordered_ind][17:59]

with PdfPages(output_pdf) as pdf:
    plt.figure(figsize=(10, 6))

    for donor, _df in df_dict.items():
        if donor in donors:
            data = _df["Entropy_mean"].dropna()
            density = gaussian_kde(data)
            x = np.linspace(data.min(), data.max(), 100)
            plt.plot(x, density(x), label=f'{donor}')

    plt.xlabel("X-axis Label")
    plt.ylabel("Density")
    plt.title("Overlayed Density Plots")
    
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small', ncol=1)
    plt.grid(True)
    
    pdf.savefig(bbox_inches='tight')
    plt.close()


output_pdf = "99_Fig/sup_fig6/002_density_cluster3.pdf"

donors = cosine_matrix.index[p.dendrogram_col.reordered_ind][59:]

with PdfPages(output_pdf) as pdf:
    plt.figure(figsize=(10, 6))

    for donor, _df in df_dict.items():
        if donor in donors:
            data = _df["Entropy_mean"].dropna()
            density = gaussian_kde(data)
            x = np.linspace(data.min(), data.max(), 100)
            plt.plot(x, density(x), label=f'{donor}')

    plt.xlabel("X-axis Label")
    plt.ylabel("Density")
    plt.title("Overlayed Density Plots")
    
    plt.legend(loc='center left', bbox_to_anchor=(1, 0.5), fontsize='small', ncol=2)
    plt.grid(True)
    
    pdf.savefig(bbox_inches='tight')
    plt.close()


# Count heterotypic connections in tiles with entropy above the cutoff
for donor, _df in df_dict.items():
    print(donor)
    print(_df[["Entropy_mean","Purity_mean"]].describe())


df_dict2 = {}

for donor, _df in df_dict.items():
    df_dict2[donor] = _df[ _df["Entropy_mean"] >= 0.5 ]


edge_df_dict2 = {}
radius = 80

for name, _df in df_dict2.items():
    points = _df[["x","y"]]
    adj_matrix = radius_neighbors_graph(X=points, radius=radius)
    edge_index = np.argwhere(adj_matrix > 0)
    edge_df = pd.DataFrame(edge_index, columns=["edge_from","edge_to"])
    edge_df["id_from"] = _df.index[edge_df["edge_from"]]
    edge_df["id_to"] = _df.index[edge_df["edge_to"]]
    edge_df["label_from"] = _df.iloc[edge_df["edge_from"], _df.columns.get_loc("label") ].values
    edge_df["label_to"] = _df.iloc[edge_df["edge_to"], _df.columns.get_loc("label")].values
    edge_df_dict2[name] = edge_df


AdjMat_dict = {}

cluster_level = ["Immune cells", "plasma", "Fibro(loose)", "Fibro(dense,regular)",
"Fibro(dense,irregular)", "lining", "vessel", "vessel(large)",
"adipose", "Stroma", "muscle", "RBC"]

for name, _edge_df in tqdm(edge_df_dict2.items()):

    neighbor_cluster_num = {}

    for cluster in cluster_level:
        edge_subset = _edge_df[ _edge_df['label_from'] == cluster]

        cluster_counts = {key: 0 for key in cluster_level}
        cluster_counts.update(edge_subset['label_to'].value_counts().to_dict())
        neighbor_cluster_num[cluster] = cluster_counts

    AdjMat = pd.DataFrame(neighbor_cluster_num)
    AdjMat = AdjMat / _edge_df.shape[0]

    AdjMat_dict[name] = AdjMat


AdjMat_dict2 = {}

for name, _edge_df in edge_df_dict2.items():

    neighbor_cluster_num = {}

    for cluster in cluster_level:
        edge_subset = _edge_df[ _edge_df['label_from'] == cluster]

        cluster_counts = {key: 0 for key in cluster_level}
        cluster_counts.update(edge_subset['label_to'].value_counts().to_dict())
        neighbor_cluster_num[cluster] = cluster_counts

    AdjMat = pd.DataFrame(neighbor_cluster_num)
    AdjMat_upper = pd.DataFrame(np.triu(AdjMat,k=1) + np.triu(AdjMat,k=1),
                                index=AdjMat.index, columns=AdjMat.columns)
    AdjMat_upper = AdjMat_upper / AdjMat_upper.sum().sum()


    AdjMat_dict2[name] = AdjMat_upper


triu_list = [ df.values[ np.triu_indices_from(df) ] for df in AdjMat_dict.values()]
pd.DataFrame(triu_list)

degree_df = pd.DataFrame(triu_list, index=AdjMat_dict.keys())
tmp = pd.DataFrame(index=cluster_level,columns=cluster_level)

for i in range(len(cluster_level)):
    for j in range(len(cluster_level)):
        term = cluster_level[i] + "_" + cluster_level[j]
        tmp.iloc[i,j] = term

columns = tmp.values[ np.triu_indices_from(tmp) ]

degree_df.columns = columns

(degree_df.sum(axis=0) == 0).value_counts()
# False    78
# Name: count, dtype: int64

degree_df.columns[np.where(degree_df.sum(axis=0) == 0)]

degree_df.loc[:,degree_df.sum(axis=0) != 0].shape
# (118, 78)

degree_df = degree_df.loc[:,degree_df.sum(axis=0) != 0]


triu_list = [ df.values[ np.triu_indices_from(df, k=1) ] for df in AdjMat_dict2.values()]
pd.DataFrame(triu_list)

degree_df2 = pd.DataFrame(triu_list, index=AdjMat_dict2.keys())
tmp = pd.DataFrame(index=cluster_level,columns=cluster_level)

for i in range(len(cluster_level)):
    for j in range(len(cluster_level)):
        term = cluster_level[i] + "_" + cluster_level[j]
        tmp.iloc[i,j] = term

columns = tmp.values[ np.triu_indices_from(tmp,k=1) ]
degree_df2.columns = columns

(degree_df2.sum(axis=0) == 0).value_counts()
# alse    66
# Name: count, dtype: int64

degree_df2.columns[np.where(degree_df2.sum(axis=0) == 0)]

degree_df2.loc[:,degree_df2.sum(axis=0) != 0].shape
# (118, 88)

degree_df2 = degree_df2.loc[:,degree_df2.sum(axis=0) != 0]


degree_df.to_csv("00_src/Degree_df_241220.csv")
degree_df2.to_csv("00_src/Degree_df2_241220.csv")

