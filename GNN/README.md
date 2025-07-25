GNN/
├── environment.yml                # conda仮想環境情報
├── README.md                      # 本ファイル
├── model_weights/
│   └── CancerCell_model_20um_5hop3sample64feat_241118.pt          # 学習済みモデルの重みファイル
└── qupath/                        # QuPathスクリプトや設定ファイルなど

# 訓練データ (SLIC tile内のラベル比率) 118 files
C:\Users\admin\Documents\東京大学 整形外科\ヒト滑膜\labels_QuPath6\labelRatio_241117
├── D001_HE.txt
├── D002_HE.txt
├── ...
└── D130_HE.txt

# 訓練データ (SLIC tile内のラベル比率) + 訓練済み特徴量のHDF5
C:\Users\admin\Documents\東京大学 整形外科\ヒト滑膜\labels_QuPath6\
└── GNN_InOut.h5

# 仮想環境の再現
```
conda env create -f environment.yml
conda activate GNN
```