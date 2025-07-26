// ----- SET THESE PARAMETERS -----
// 書き出し対象のアノテーションクラス
List classNames = ["Immune cells", "plasma", "Fibro(loose)", "Fibro(dense,regular)",
"Fibro(dense,irregular)", "lining", "vessel", "vessel(large)",
"adipose", "Stroma", "muscle", "RBC"]

int downsample = 1

// タイル中に含まれるアノテーション領域の比。labelAreaRatioより少ないアノテーション面積であればそのタイル画像はかき出さない。
double labelAreaRatio = 0
int tileSize = 512  // タイルの辺の長さ px
int tileOverlap = 256  // タイルの重複幅 px
def imageExtension = ".jpg"
def pathOutput = "C:/Users/admin/Documents/東京大学 整形外科/ヒト滑膜/hSyn_QuPath3//Train6_241107"
// --------------------------------

import ij.*
import ij.IJ
import ij.measure.ResultsTable

def imageData = getCurrentImageData()
def server = imageData.getServer()
server.setMetadata(server.getOriginalMetadata())
def imageName = GeneralTools.getNameWithoutExtension(server.getMetadata().getName())

print "===== ${imageName} ====="

// 保存先フォルダの作成
def imageOutDir = buildFilePath(pathOutput, "Images")
def labelOutDir = buildFilePath(pathOutput, "Labels")
mkdirs(pathOutput); mkdirs(imageOutDir); mkdirs(labelOutDir)

// annotaionのLabel serverを構築
def tempServer = new LabeledImageServer.Builder(imageData)
    .backgroundLabel(0, ColorTools.WHITE)
    .downsample(downsample) 
    .multichannelOutput(false)
    .grayscale()

// Label serverにアノテーションクラスを登録
def counter = 1
classNames.each { currClassName ->
    tempServer.addLabel(currClassName, counter)  // Choose output labels (the order matters!)
    counter++;
}
// finally, build server
def labelServer = tempServer.build()
print "Label server contains: ${labelServer.getLabels()}"

// タイル辺長と重複辺長をダウンサンプリング係数で補正
tileSize = tileSize*downsample
tileOverlap = tileOverlap*downsample

// 選択画像のserver情報から画像pxサイズを取得
int width = getCurrentServer().getWidth(); print "Image px width: ${width}"
int height = getCurrentServer().getHeight(); print "Image px height: ${height}"

def plane = ImagePlane.getDefaultPlane()

// タイル辺長と重複辺長から四角ROIのROIクラスを作っていく。
def coordList = []
def nextStart = tileSize - tileOverlap // 次のタイル開始位置情報に変換

for(int i=0; i<width; i+=nextStart){
    for(int j=0; j<height; j+=nextStart){
        coordList << [i,j]
    }
}

tileList = coordList.collect{ROIs.createRectangleROI(it[0], it[1], tileSize, tileSize, plane)}

print "${tileList.size()} tiles can be created"
def serverPath = server.getPath()
def labelServerPath = labelServer.getPath()

// タイルを1つずつ処理
def currNum = 0
tileList.each{
    // アノテーションのserverからタイルROIの場所をリクエストする
    def requestLabelROI = RegionRequest.createInstance(labelServerPath, downsample, it)
    // RegionRequestクラスからBufferedImageクラスを生成
    def label = labelServer.readRegion(requestLabelROI) // class java.awt.image.BufferedImage
    
    // BufferedImageからImageJのImagePlusクラスに変換
    label = new ImagePlus("Mask", label) // class ij.ImagePlus
    
    // ImageJのMeasureでアノテーション面積を測定
    ij.IJ.run("Set Measurements...", "area_fraction redirect=None decimal=3");
    ij.IJ.run(label, "Measure", "");
    def maskArea = ResultsTable.getResultsTable().getValue("%Area",0)
    
    // マスク画像のマスクエリア割合が閾値より大きければ書き出し
    if (maskArea >= labelAreaRatio*100) {
        // 書き出しファイル名に入れる座標情報作成
        int h = it.getBoundsHeight(); int w = it.getBoundsWidth()
        int x = it.getBoundsX(); int y = it.getBoundsY()
        def tileInfo = "[d=${downsample},x=${x},y=${y},w=${w},h=${h}]"   
        
        // 書き出しパス作成
        def saveLabelName = "${imageName} ${tileInfo}.png"        
        def saveLabelPath = buildFilePath(labelOutDir, saveLabelName)
        // RegionRequestクラスを使ってアノテーションのserverからマスク画像を書き出し
        writeImageRegion(labelServer, requestLabelROI, saveLabelPath)        
        
        // 明視野画像の方も書き出し
        def requestROI = RegionRequest.createInstance(serverPath, downsample, it)        
        def saveName = "${imageName} ${tileInfo}${imageExtension}"
        def savePath = buildFilePath(imageOutDir, saveName)
        writeImageRegion(server, requestROI, savePath)        
        
        currNum += 1
    }    
}
print "${currNum} tiles were passed threshold!!"
