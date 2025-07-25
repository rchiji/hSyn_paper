// Parameters ----------------------------------------------------
// ファイルへのパス
String path = buildFilePath(PROJECT_BASE_DIR, "SLICTile_cluster_241126",
    getCurrentImageName().split("_HE")[0] + ".txt")
// 区切り文字
String delim = "\t"
// 対象Object一覧
def objects = getTileObjects()
// ---------------------------------------------------------------

import qupath.lib.gui.measure.ObservableMeasurementTableData

// 1行ずつ読み込み
def lines = new File(path).readLines()

// ヘッダー行を取り出し
def header = lines.pop().split(delim)

// Object IDとObjectのMapを作成
Map objectsById = objects.groupBy(object -> object.getID().toString())

// 1行ずつ処理
for (def line in lines) {
    // Mapに項目名と値を登録
    Map map = [:]
    def content = line.split(delim) // 区切り文字で要素を分ける
    for (int i = 0; i < header.size(); i++) {
        map[header[i]] = content[i]
        }

    // Object IDを元に対象Objectを探す
    String id = map['ID']
    def object = objectsById[id][0]
    
    if (object == null) {
        continue
    }
    def classification = getPathClass("KMeans: " + map["Cluster"])
    object.setPathClass(classification)
    }
