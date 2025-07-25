// Parameter ---------------
double downsample = 5 // ダウンサンプリング係数
String extension = ".ome.tiff" // 保存時の拡張子
String saveDir = "label_ometiff_model7_241113" // 保存先のフォルダ
// -------------------------

import qupath.lib.images.servers.LabeledImageServer

def imageData = getCurrentImageData()

// 保存先のフォルダを作成
String imageName = getCurrentImageNameWithoutExtension()
String pathOutput = buildFilePath(PROJECT_BASE_DIR, saveDir)
mkdirs(pathOutput)

// PathClass一覧を取得
List classNames = [
    "Immune cells", 
    "plasma",
    "Fibro(loose)",
    "Fibro(dense,regular)",
    "Fibro(dense,irregular)",
    "lining", 
    "vessel",
    "vessel(large)",
    "adipose",
    "Stroma",
    "muscle",
    "RBC"
    ]

// LabeledImageServerを構築
def labelServer = new LabeledImageServer.Builder(imageData)
for(i = 0; i < classNames.size(); i++) {
    labelServer.addLabel(classNames[i], i+1)
    }
labelServer = labelServer
    .downsample(downsample)
    .multichannelOutput(false)
    .build()

// 書き出しエリア　（Z stack画像があるため、Z=0の領域を明示的に指定）
def request = RegionRequest.createInstance(labelServer.getPath(), downsample,
                                            0,0, labelServer.getWidth(), labelServer.getHeight(), 0,0)

// 書き出し
String savePath = "${pathOutput}/${imageName}${extension}"
writeImageRegion(labelServer, request, savePath) 

