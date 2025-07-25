double largeObjectArea = 10000 // Area µm^2

print getCurrentImageName()


def objects = getAnnotationObjects()

// ----- 1. 削除判定する対象Objectの取得 -----
def checkObjects = objects.findAll {
    it.getPathClass().toString() == "adipose"
}

print "Check objects num: " + checkObjects.size()


// ----- 2. 距離計測に使用する大きめのObjectの取得 -----
import qupath.lib.gui.measure.ObservableMeasurementTableData

// ObservableMeasurementTableDataインスタンスを作成
def ob = new ObservableMeasurementTableData()

// obにアノテーションを登録
ob.setImageData(getCurrentImageData(), objects)

def largeObjects = objects.findAll {
    ob.getNumericValue(it, "Area µm^2") > largeObjectArea
}


// ----- 3. Object間の距離でフィルタリング -----

// Object間の距離の記録先Map
Map distMap = new HashMap( checkObjects.size() )
checkObjects.each {
    distMap[it] = new HashMap( largeObjects.size() )
}


// 総当たりで距離計算を並列処理
checkObjects.parallelStream().forEach { object1 ->
    def ROI1 = object1.getROI()
    largeObjects.parallelStream().forEach { object2 ->
        if ( object1 != object2 ) {  // 同じオブジェクト同士の計算は除外
            double distance = RoiTools.getBoundaryDistance( ROI1, object2.getROI() )
            synchronized (distMap) {  // distMapへのアクセスは同期させる
                distMap[object1][object2] = distance
                }
            }
        }
    }

// 大きいObjectと接していないAdipoは削除
def deleteObjects = checkObjects.parallelStream().filter {
    distMap[ it ].values().min() > 0
    }.collect{it}

print "Delete objects num: " + deleteObjects.size()

removeObjects( deleteObjects, true )