# Project Structure and Data Organization
```
GNN/
├── environment.yml                                          # Conda environment file
├── README.md                                                # Project documentation
├── 2401118_SLICTile_GNN.ipynb                               # Training record notebook
├── model_weights/                                           # Trained model weights
│   └── CancerCell_model_20um_5hop3sample64feat_241118.pt
|
├── qupath/                                                  # QuPath project for processing WSI prediction label images (imported from label_ometiff_model7_241113)
|   └── scripts/                                             # Scripts for generating SLIC superpixels and mapping GNN features back to QuPath
|
└── data/
    ├── labelRatio_241117/                                   # Label composition per SLIC tile (download required)
    |   ├── D001_HE.txt
    |   ├── ...
    |   └── D130_HE.txt
    |
    ├── SLICTile_Features_241125/                            # GNN-derived features per SLIC tile (download required)
    │   ├── D001.txt
    │   ├── ...
    │   └── D130.txt
    │
    └── GNN_InOut.h5                                         # GNN input/output data (labels + features, HDF5 format, download required)

```

# Data Availability
The datasets used in this study — including label composition within SLIC tiles, 64-dimensional SLIC tile features, and GNN-derived input/output data (labels + features in HDF5 format) — are available via Zenodo:
https://zenodo.org/communities/human-synovial-histology

### Contents
- **Label composition files (per SLIC tile, 118 files)**  
 -> *SLIC tile label composition features for human synovial histology 01-02* 

- **GNN-derived features (per SLIC tile, 118 files)**  
 -> *Graph neural network–derived SLIC tile features for human synovial histology 01-02*

- **GNN input/output data (labels + features, HDF5 format)**  
 -> *Graph neural network input and output data for human synovial histology analysis*


# Environment Setup
```
conda env create -f environment.yml
conda activate GNN
```
