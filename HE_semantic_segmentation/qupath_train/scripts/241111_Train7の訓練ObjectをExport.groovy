String saveDir = buildFilePath(PROJECT_BASE_DIR, 'Train7_trainObjects_241111')
mkdirs(saveDir)

// geojsonファイルを保存
String saveName = getCurrentImageNameWithoutExtension() 
String savePath = buildFilePath(saveDir, "${saveName}.zip")
exportAllObjectsToGeoJson(savePath)