// Parameters ------
// ダウンサンプリング係数
double downsample = 1
// 対象Object一覧
def objects = getTileObjects()

String saveDir = "labelRatio_241117"
mkdirs buildFilePath(PROJECT_BASE_DIR, saveDir)

// 計測結果の書き出し先
String imageName = getCurrentImageNameWithoutExtension()
String filePath = buildFilePath(PROJECT_BASE_DIR, saveDir, imageName + ".txt")

//
List className = ["Background","Immune cells", "plasma", "Fibro(loose)", "Fibro(dense,regular)",
"Fibro(dense,irregular)", "lining", "vessel", "vessel(large)",
"adipose", "Stroma", "muscle", "RBC"]
// -----------------

import ij.ImagePlus
import ij.IJ

def server = getCurrentServer()

List labelTables = []
objects.each {
    // 対象領域情報を作成
    def roi = it.getROI()
    
    // ROIの左上端の座標
    int xPosition = roi.getBoundsX()
    int yPosition = roi.getBoundsY()
    
    // 解析対象情報作成
    def request = RegionRequest.createInstance(server.getPath(), downsample, roi)
    
    // ImagePlusの画像を作成
    ImagePlus imp = IJTools.convertToImagePlus(server, request).getImage()
    
    int width = imp.getWidth()
    int height = imp.getHeight()

    // QuPath ROIをImageJ ROIに変換
    def roiIJ = IJTools.convertToIJRoi(roi,-roi.getBoundsX(),-roi.getBoundsY(),1)
    
    // ImagePlusにRoiを登録
    imp.setRoi(roiIJ)
    
    // 結果の集計先Map
    Map labelTable = [:]
    labelTable["ID"] = it.getID()
    labelTable["x"] = roi.getCentroidX()
    labelTable["y"] = roi.getCentroidY()
    
    // 背景以外のラベル番号の初期値を0に設定
    (1..12).each { label ->
        labelTable.put(label,0)
        }
    
    // 座標は0始まりなので0からPx数-1の範囲を集計
    for( int x in 0..width-1 ) {
       for( int y in 0..height-1 ) {
           // XY座標がROIの範囲内の時のみ集計
           if ( roiIJ.contains(x,y) ) {
               label = imp.getPixel(x,y)
               // 背景輝度でなければ集計（PNG画像が4ch画像として扱われており[8,0,0,0]のように0ch目がラベル、2-4ch目は0のリストが返ってくる。）
               if (label[0] > 0) {
                   // 集計先に1を加算
                   labelTable[label[0]] += 1               
                   } // <- if (label[0] > 0)               
               } // <- if ( roiIJ.contains(x,y) )
           } // <- for( int y in 0..height-1 )
       } // <- for( int x in 0..width-1 )
    
    labelTables << labelTable
    
    // ID, X, Y列を除いたリスト取得
    List<Number> values = labelTable.values().toArray()[3..(labelTable.size()-1)]
    int maxValue = values.max()
    int key = labelTable.find{ key,value ->
        value == maxValue
        }.key

    it.setName(className[key])
    def ml = it.getMeasurementList()
    ml.put("Max label: ",key)
    
    values.eachWithIndex{ v, i ->
        ml.put("Px num: " + className[i+1], v)
        }
    }
    
    
print labelTables

// ----- 書き出し -----
print "Writing..."

new File(filePath).withWriter { writer ->
    def bufferedWriter = new BufferedWriter(writer)
    keys = labelTables[0].keySet()
    for(String key:keys) {
        bufferedWriter.write(key + "\t")
        }
    bufferedWriter.write("\n")
    
    StringBuilder sb = new StringBuilder()
    labelTables.each { labelTable ->
        Set<String> keys = labelTable.keySet()
        
        sb.setLength(0)  // Clear the StringBuilder
        for(key:keys) {
            sb.append(labelTable.get(key))
            sb.append("\t")
            }
        sb.append('\n')
        bufferedWriter.write(sb.toString())
    }
    bufferedWriter.flush()
}

print "Done!!"    