// Parameters ------
// 対象Object一覧
def objects = getTileObjects()

String saveDir = "tile_label_coord_241115"
mkdirs buildFilePath(PROJECT_BASE_DIR, saveDir)

// 計測結果の書き出し先
String imageName = getCurrentImageNameWithoutExtension()
String filePath = buildFilePath(PROJECT_BASE_DIR, saveDir, imageName + ".txt")
// -----------------

List labelTables = []
objects.each {
    
    // 対象領域情報を作成
    def roi = it.getROI()
        
    Map labelTable = [:]
    labelTable["ID"] = it.getID()
    labelTable["x"] = roi.getCentroidX()
    labelTable["y"] = roi.getCentroidY()
    labelTable["label"] = it.getPathClass().toString()

    
    labelTables << labelTable
    }
    
print labelTables  
// ----- 書き出し -----
print "Writing..."

new File(filePath).withWriter { writer ->
    def bufferedWriter = new BufferedWriter(writer)
    def keys = labelTables[0].keySet()
    
    // ヘッダーの書き込み
    keys.eachWithIndex { key, idx ->
        bufferedWriter.write(key)
        if (idx < keys.size() - 1) bufferedWriter.write("\t")  // 最後の列にタブを追加しない
    }
    bufferedWriter.write("\n")
    
    // データの書き込み
    StringBuilder sb = new StringBuilder()
    labelTables.each { labelTable ->
        sb.setLength(0)  // Clear the StringBuilder
        keys.eachWithIndex { key, idx ->
            sb.append(labelTable.get(key))
            if (idx < keys.size() - 1) sb.append("\t")  // 最後の列にタブを追加しない
        }
        sb.append("\n")
        bufferedWriter.write(sb.toString())
    }
    bufferedWriter.flush()
}


print "Done!!"    
