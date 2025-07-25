import qupath.lib.gui.measure.ObservableMeasurementTableData
// import qupath.lib.gui.commands.SummaryMeasurementTableCommand

// 書き出し対象のObject
List objects = getAnnotationObjects()
// MeasurementListの項目一覧を取得

// ObservableMeasurementTableDataのインスタンス作成
def ob = new ObservableMeasurementTableData()
// 書き出し対象のObjectを登録
ob.setImageData(getCurrentImageData(), objects)
List terms = ob.getAllNames()

// 保存ディレクトリのパス
String saveDir = buildFilePath(PROJECT_BASE_DIR, "Annotation_ratio_241204")
mkdirs saveDir
// 計測結果の書き出し先
String imageName = getCurrentImageNameWithoutExtension()
String filePath = buildFilePath(saveDir, imageName + ".txt")

// Fileクラスの変数作成
File file = new File(filePath)

// 同じファイルがあった場合に削除
file.delete()  // この行を削除すると追記になる
// ファイル作成
file.createNewFile()


// StringBuilderを利用して内容をバッファに格納
StringBuilder buffer = new StringBuilder()

// header部分を書いておく
buffer.append("Object ID" +"\t")
terms.eachWithIndex { key, i ->
    buffer.append(key)
    if ( i < terms.size() -1 ) {
        buffer.append("\t") 
        } else {
        buffer.append("\n")
        }
    }

objects.each {
    // 書き出すObject IDを取得
    String id = it.getID()
    buffer.append( "${id}\t" ) // 追記
    
    // 書き出す内容を取得
    def ml = it.getMeasurementList()
    
    // 計測項目を取り出して追記
    terms.eachWithIndex { key, i ->
        def value = ob.getNumericValue(it, key)
        if ( value.isNaN() ) {
            value = ob.getStringValue(it, key)
            }
        
        buffer.append( value )
        if ( i < terms.size() -1 ) { 
            buffer.append("\t")
            } else {
            buffer.append("\n")
            }
        }
    }



file.withWriter { writer ->
    writer.write( buffer.toString() )
    }