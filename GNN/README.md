# 📦 プロジェクト構成とデータ配置
```
GNN/
├── environment.yml                                  # conda仮想環境ファイル
├── README.md                                        # 本ファイル
├── model_weights/                                   # 学習済みモデル
│   └── CancerCell_model_20um_5hop3sample64feat_241118.pt
└── qupath/                                          # QuPath用のスクリプトや設定ファイル

├── data/
│   ├── labelRatio/                           # SLICタイルごとのラベル比率ファイル（118個） <-- need download
│   │   ├── D001_HE.txt
│   │   ├── ...
│   │   └── D130_HE.txt

│   └── GNN_InOut.h5                                 # 訓練済み特徴量 + ラベル情報（HDF5形式） <-- need download

```

# 📦 ダウンロードリンク一覧

- **ラベル比率ファイル（SLIC tile単位, 118 files）**  
  [labelRatio フォルダ](https://drive.google.com/drive/folders/1G-Y3E7iS7MdGGhgSz1ogzkP2WI6Z1wm0?usp=drive_link)

- **訓練済み特徴量（ラベル＋特徴量, HDF5形式）**  
  [GNN_InOut.h5](https://drive.google.com/file/d/1riQMyqboWnH4b2S0tbYj88Ww5YaTtco9/view?usp=drive_link)


# 📦 仮想環境の再現
```
conda env create -f environment.yml
conda activate GNN
```