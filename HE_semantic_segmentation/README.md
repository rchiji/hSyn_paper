HE_semantic_segmentation/
├── environment.yml                # conda仮想環境情報
├── README.md                      # 本ファイル
├── model_weights/                 # 学習済みモデルの重みファイル
│   └── model7_241111.h5
└── qupath/                        # QuPathスクリプトや設定ファイルなど

## model_weight
https://drive.google.com/file/d/19RXABNNg3Ww48RUgZqlIS-PtlxSF2HeV/view?usp=drive_link

# HE画像ファイルパス
C:\Users\admin\Documents\東京大学 整形外科\ヒト滑膜\HE\
├── D001_HE.ndpi
├── D002_HE.ndpi
├── ...
└── D130_HE.ndpi

# 訓練データ
C:\Users\admin\Documents\東京大学 整形外科\ヒト滑膜\hSyn_QuPath3\Train7_241111
├── images/ # 2690 files
|   ├── D001_HE [d=5,x=49282,y=27673,w=2899,h=1652] [d=1,x=0,y=0,w=512,h=512].jpg
|   ├── ...
│   └── D130_HE [d=5,x=9661,y=17238,w=4628,h=4764] [d=1,x=768,y=768,w=512,h=512].jpg
└── labels/ # 2690 files
|   ├── D001_HE [d=5,x=49282,y=27673,w=2899,h=1652] [d=1,x=0,y=0,w=512,h=512].png
|   ├── ...
│   └── D130_HE [d=5,x=9661,y=17238,w=4628,h=4764] [d=1,x=768,y=768,w=512,h=512].png

# 予測後データ
C:\Users\admin\Documents\東京大学 整形外科\ヒト滑膜\hSyn_QuPath3\
├── label_ometiff_model7_241113\                 # ラベルマップ（OME-TIFF形式）
├── Prediction_label_json_model7_241113\        # 推論結果のJSON形式
└── rendered_thumnails_model7_241113\           # サムネイル可視化画像

# 仮想環境の再現
```
conda env create -f environment.yml
conda activate tf
```