// ----- Paramters ----------------
// この面積より小さいObjectは不採用
double minAreaThreshold = 5000

// この面積より大きいObjectは必ず採用
double largeObjectArea = 100000

// px。minArea閾値を超えても、距離閾値以内にlargeObjectAreaを満たすObjectが無いと不採用
double distThreshold = 1000 
// ---------------------------------

selectAnnotations();
runPlugin('qupath.lib.plugins.objects.SplitAnnotationsPlugin', '{}')


print getCurrentImageName()

// ----- 1. 対象Objectの取得 -----
def objects = getAnnotationObjects().findAll {
    it.getPathClass().toString() == "Region*"
}
removeObjects( objects, true)

// ----- 2. Objectの面積閾値でフィルタリング ----
import qupath.lib.gui.measure.ObservableMeasurementTableData

// ObservableMeasurementTableDataインスタンスを作成
def ob = new ObservableMeasurementTableData()

// obにアノテーションを登録
ob.setImageData(getCurrentImageData(), objects)

print "Befor filtering: " + objects.size()
objects = objects.findAll {
    ob.getNumericValue(it, "Area µm^2") > minAreaThreshold
}
print "After size filtering: " + objects.size()


// ----- 3. Object間の距離でフィルタリング -----

// Object間の距離の記録先Map
Map distMap = new HashMap(objects.size())
objects.each {
    distMap[it] = new HashMap(objects.size())
}

Map roiMap = objects.collectEntries { [(it): it.getROI()] }
Map areaMap = objects.collectEntries { [(it): ob.getNumericValue(it, "Area µm^2") ] }

// 総当たりで距離計算を並列処理
objects.parallelStream().forEach { object1 ->
    objects.parallelStream().forEach { object2 ->
        if (object1 != object2) {  // 同じオブジェクト同士の計算は除外
            double distance = RoiTools.getBoundaryDistance(roiMap[object1], roiMap[object2])
            synchronized (distMap) {  // distMapへのアクセスは同期させる
                distMap[object1][object2] = distance
                distMap[object2][object1] = distance
                }
            }
        }
    }


// distThresholdの距離に閾値以上の面積を持つRegion*が無いものは除く （周囲が小さいものしかないもの）
objects = objects.parallelStream().filter {
    
    // Objectのサイズが大きかったら周囲関係なしに採用
    if ( areaMap[it] > largeObjectArea ) {
        return true
        }
        
    // 最も近いObjectが1 mm以上離れていたら不採用
    else if ( distMap[ it ].values().min() > distThreshold ) {
        return false
        }
        
    // 1 mm以内にあるObjectが小さいものばかりなら不採用
    else {
        // 1000 px以内の距離のObjectのみに限定
        Map distMap2 = distMap[ it ].findAll{ key,value -> value <= distThreshold }
        // 近傍Objectの面積を取得
        List nearObjectAreas = distMap2.keySet().collect{ areaMap[it] }
        // 近傍Objectの面積が大きいか否かを返す
        return nearObjectAreas.max() > largeObjectArea
        }
    
    }.collect{it}

print "After distance filtering: " + objects.size()


def finalRegionROI = RoiTools.union( objects.collect{ roiMap[it] } )
def hierarchy = getCurrentHierarchy()
def remainObjects = hierarchy.getObjectsForROI( objects[0].class, finalRegionROI )

clearAnnotations()
addObjects( remainObjects )