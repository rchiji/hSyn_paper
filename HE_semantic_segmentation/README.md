# 📦 プロジェクト構成とデータ配置
```
HE_semantic_segmentation/
├── environment.yml                         # conda仮想環境ファイル
├── README.md                               # 説明ドキュメント
├── model_weights/                          # 学習済みモデルの重みファイル
│   └── model7_241111.h5                    # <-- need download
├── qupath/                                 # QuPath用のスクリプトや設定ファイル

├── data/
│   ├── HE/                             # HEスライド画像（NDPI)
│   │   ├── D001_HE.ndpi
│   │   ├── ...
│   │   └── D130_HE.ndpi

│   ├── Train/                      # 訓練用パッチ画像およびマスク <-- need download
│   │   ├── images/ # 2690 files
│   │   │   ├── D001_HE [d=5,x=49282,y=27673,w=2899,h=1652] [d=1,x=0,y=0,w=512,h=512].jpg
│   │   │   ├── ...
│   │   │   └── D130_HE [d=5,x=9661,y=17238,w=4628,h=4764] [d=1,x=768,y=768,w=512,h=512].jpg
│   │   └── labels/ # 2690 files
│   │       ├── D001_HE [d=5,x=49282,y=27673,w=2899,h=1652] [d=1,x=0,y=0,w=512,h=512].png
│   │       ├── ...
│   │       └── D130_HE [d=5,x=9661,y=17238,w=4628,h=4764] [d=1,x=768,y=768,w=512,h=512].png
│   │

│   └── predictions/                         # モデルの推論結果群
│       ├── label_ometiff_model7_241113/     # OME-TIFF形式ラベル画像
│       ├── Prediction_label_json_model7_241113/  # JSON形式ラベルデータ
│       └── rendered_thumnails_model7_241113/ # 可視化サムネイル

```

# 📦 ダウンロードリンク一覧

- **モデル重み**  
  [model7_241111.h5](https://drive.google.com/file/d/19RXABNNg3Ww48RUgZqlIS-PtlxSF2HeV/view?usp=drive_link)

- **HE画像（NDPIスライド）**  
  [HEフォルダ](https://drive.google.com/drive/folders/1zCC_lSDlx4vzBxic3E_GgF6Djyitz-ZB?usp=drive_link)

- **訓練データ（パッチ画像・マスク）**  
  [Train7_241111](https://drive.google.com/drive/folders/18ybfXBofi4t-X2OJG2F0aTgNp9FRiOOZ?usp=drive_link)

- **推論結果データ**  
  - [label_ometiff_model7_241113](https://drive.google.com/drive/folders/19FzmIiERxySgiVUBLJM1ePCDAJIFRlu1?usp=drive_link)  
  - [Prediction_label_json_model7_241113](https://drive.google.com/drive/folders/19EOG1G_A-av6Y4qssU5p1B4vcmk2epr6?usp=drive_link)  
  - [rendered_thumnails_model7_241113](https://drive.google.com/drive/folders/197165woG291TywZgyKePgV1H16TPVSjE?usp=drive_link)


# 仮想環境の再現
```
conda env create -f environment.yml
conda activate tf
```