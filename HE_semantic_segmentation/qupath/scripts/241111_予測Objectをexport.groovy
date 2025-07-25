String saveDir = buildFilePath(PROJECT_BASE_DIR, 'Prediction_label_json_241018')
mkdirs(saveDir)

// geojsonファイルを保存
String saveName = getCurrentImageNameWithoutExtension() 
String savePath = buildFilePath(saveDir, "${saveName}.zip")
exportAllObjectsToGeoJson(savePath)