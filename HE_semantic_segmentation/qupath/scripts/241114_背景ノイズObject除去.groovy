// --- 1. annotation objectを分割 ---
selectAnnotations();
runPlugin('qupath.lib.plugins.objects.SplitAnnotationsPlugin','{}')
resetSelection()

// --- 2. 背景にIgnore* object作成 ---

// Parameter -------------------------------------------
// ① ImageDataOp
def op = ImageOps.buildImageDataOp()
                 .appendOps(
                     ImageOps.Channels.mean(),
                     ImageOps.Filters.gaussianBlur(5),
                     ImageOps.Threshold.threshold(220)
                     )
                

// ② ダウンサンプリング係数
double downsample = 64

// ③ 二値化の輝度値とPathClassの対応表
Map classifications = Map.of(
    1, getPathClass("Ignore*")
    )

// Object作成時のオプション
double minArea = 10000000
double minHoleArea = 0
boolean INCLUDE_IGNORED = true
boolean DELETE_EXISTING = false
boolean SELECT_NEW = false
boolean SPLIT = false
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


// --- 3. Ignore*に内包されるObjectを削除 ---
def ignoreROI = getAnnotationObjects().find {
    it.getPathClass().toString() == "Ignore*"
    }.getROI()

def hierarchy = getCurrentHierarchy()
List deleteObjects = hierarchy.getObjectsForROI( null, ignoreROI )

removeObjects(deleteObjects, true)

// Z軸があるものは各stackに作られているものが残るので追加削除
def ignores = getAnnotationObjects().findAll {
    it.getPathClass().toString() == "Ignore*"
    }
if( ignores.size() > 0) {
    removeObjects( ignores, true)
    }
