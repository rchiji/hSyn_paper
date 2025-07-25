```
GNN/
├── environment.yml                # conda仮想環境情報
├── README.md                      # 本ファイル
├── model_weights/
│   └── CancerCell_model_20um_5hop3sample64feat_241118.pt          # 学習済みモデルの重みファイル
└── qupath/                        # QuPathスクリプトや設定ファイルなど
```

# 訓練データ (SLIC tile内のラベル比率) 118 files
```
C:\Users\admin\Documents\東京大学 整形外科\ヒト滑膜\labels_QuPath6\labelRatio_241117
├── D001_HE.txt
├── D002_HE.txt
├── ...
└── D130_HE.txt
```
https://drive.google.com/drive/folders/1G-Y3E7iS7MdGGhgSz1ogzkP2WI6Z1wm0?usp=drive_link

# 訓練データ (SLIC tile内のラベル比率) + 訓練済み特徴量のHDF5
```
C:\Users\admin\Documents\東京大学 整形外科\ヒト滑膜\labels_QuPath6\
└── GNN_InOut.h5
```
https://drive.google.com/file/d/1riQMyqboWnH4b2S0tbYj88Ww5YaTtco9/view?usp=drive_link

# 仮想環境の再現
```
conda env create -f environment.yml
conda activate GNN
```