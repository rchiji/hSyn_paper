# hSynovium image analysis
## Overview
This repository provides code and data for the analysis of synovial tissue architecture in knee osteoarthritis (OA) using a combination of deep learning–based semantic segmentation and graph neural network (GNN) modeling.
We developed a computational pathology framework to analyze hematoxylin and eosin (H&E)-stained whole-slide images (WSIs) of human synovium. The pipeline integrates:
1. Semantic segmentation of histological structures  
2. Spatial analysis of tissue organization  
3. Graph neural network–based modeling of microenvironmental interactions  

---

## Module Description

### hSynovium H&E WSI semantic segmentation
Directory: `HE_semantic_segmentation/`

- Performs semantic segmentation of H&E-stained WSIs  
- Extracts histological structures such as:
  - Synovial lining
  - Immune cell infiltration
  - Microvessels
  - Fibrotic regions  
- Outputs spatial label maps used for downstream analysis

---

### GNN
Directory: `GNN/`

- Constructs SLIC-based graphs from segmentation outputs  
- Models spatial relationships between tissue components  

---

### Postprocess
Directory: `Postprocess/`

- Performs downstream analysis based on semantic segmentation and GNN outputs  
- Computes quantitative metrics and spatial features  
- Generates figures used in the manuscript 