 // Parameter -------------------------------------------
// ① ImageDataOp
def op = ImageOps.buildImageDataOp()
                 .appendOps(
                     ImageOps.Threshold.threshold(0)
                     )
                

// ② ダウンサンプリング係数
double downsample = 1

// ③ 二値化の輝度値とPathClassの対応表
Map classifications = Map.of(
    1, getPathClass("Region*")
    )

// Object作成時のオプション
double minArea = 0
double minHoleArea = 0
boolean INCLUDE_IGNORED = true
boolean DELETE_EXISTING = false
boolean SELECT_NEW = false
boolean SPLIT = true
// -----------------------------------------------------

import qupath.opencv.ml.pixel.PixelClassifiers

// ImageData
def imageData = getCurrentImageData()

// 解像度
def resolution = imageData.getServer().getPixelCalibration().createScaledInstance(downsample,downsample)

// ピクセル分類条件を作成
def classifier = PixelClassifiers.createClassifier(op, resolution, classifications)


// オプションをparameterの指定に応じて変更させながらリスト化
List options = [
    INCLUDE_IGNORED == true ? PixelClassifierTools.CreateObjectOptions.INCLUDE_IGNORED : null,
    DELETE_EXISTING == true ? PixelClassifierTools.CreateObjectOptions.DELETE_EXISTING : null,
    SELECT_NEW == true ? PixelClassifierTools.CreateObjectOptions.SELECT_NEW : null,
    SPLIT == true ? PixelClassifierTools.CreateObjectOptions.SPLIT : null    
    ]

// ピクセル分類の実行
PixelClassifierTools.createAnnotationsFromPixelClassifier(
    imageData,
    classifier,
    minArea,
    minHoleArea,
    *options)