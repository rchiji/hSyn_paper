String saveDir = buildFilePath(PROJECT_BASE_DIR, 'Train6_label_241107')
mkdirs(saveDir)

// geojsonファイルを保存
String saveName = getCurrentImageNameWithoutExtension() 
String savePath = buildFilePath(saveDir, "${saveName}.zip")
exportAllObjectsToGeoJson(savePath)