# Project Structure and Data Organization
```
HE_semantic_segmentation/
├── environment.yml                                   # Conda environment file
├── README.md                                         # Project documentation
├── 241111_Train7_model_training.ipynb                # Training record notebook
├── Custom_semantic_keras/                            # Custom scripts used for training
├── model_weights/                                    # Trained model weights
|   └── model7_241111.h5                              # Download required (see below)
|
├── qupath/                                           # QuPath project for H&E image processing
|   ├── predictions/                                  # Model prediction results (download required)
|   |       ├── label_ometiff_model7_241113/          # OME-TIFF label images
|   |       ├── Prediction_label_json_model7_241113/  # JSON label data
|   |       └── rendered_thumbnails_model7_241113/    # Visualization thumbnails
|   |
|   └── scripts/                                      # Scripts for importing predictions into QuPath and post-processing
|
├── qupath_train/                                     # QuPath project for annotation
|   ├── train_sources/                                # Source images for annotation
|   └── Train7_trainObjects_241111/                   # Ground truth labels (GeoJSON)
|
└── data/
    ├── HE/                                           # H&E whole slide images (NDPI) (download required)
    |   ├── D001_HE.ndpi
    |   ├── ...
    |   └── D130_HE.ndpi
    |
    └── Train/                                        # Training patches and masks (download required)
        ├── images/ # 2690 files
        |   ├── D001_HE [d=5,x=49282,y=27673,w=2899,h=1652] [d=1,x=0,y=0,w=512,h=512].jpg
        |   ├── ...
        |   └── D130_HE [d=5,x=9661,y=17238,w=4628,h=4764] [d=1,x=768,y=768,w=512,h=512].jpg
        └── labels/ # 2690 files
            ├── D001_HE [d=5,x=49282,y=27673,w=2899,h=1652] [d=1,x=0,y=0,w=512,h=512].png
            ├── ...
            └── D130_HE [d=5,x=9661,y=17238,w=4628,h=4764] [d=1,x=768,y=768,w=512,h=512].png

```


# Data Availability
The datasets used in this study — including H&E whole slide images (NDPI), training data (patch images and masks), trained model weights, and prediction results — are available via Zenodo:
https://zenodo.org/communities/human-synovial-histology

### Contents
- **H&E whole slide images (NDPI)**  
 -> *Human synovial H&E whole slide images 01–04*

- **Training data (patch images and masks)**  
 -> *Training patch images and masks for semantic segmentation of human synovial H&E whole slide images*

- **Model weights**  
 -> *Deep learning model for semantic segmentation of human synovial histology*

- **Prediction results**  
 -> *Rendered segmentation thumbnails for human synovial H&E whole slide images 01–02*
 -> *Semantic segmentation objects (GeoJSON) for human synovial H&E whole slide images 01–02*
 -> *Semantic segmentation labels (OME-TIFF) for human synovial H&E whole slide images 01–02*


# Environment Setup
```
conda env create -f environment.yml
conda activate tf
```


# Training
Training was performed using scripts in the Custom_semantic_keras directory.
For details, refer to: 241111_Train7_model_training.ipynb


# Inference
Inference is performed by executing Python scripts from within QuPath.
Refer to: qupath/scripts/01_semantic_segmentation_full_slide_model7.groovy 
The directory `qupath/scripts/custom_scripts/` contains Python packages used within QuPath.