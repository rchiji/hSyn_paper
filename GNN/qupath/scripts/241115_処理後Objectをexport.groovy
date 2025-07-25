String saveDir = buildFilePath(PROJECT_BASE_DIR, 'PostProcessed_objects_json_241116')
mkdirs(saveDir)

// geojsonファイルを保存
String saveName = getCurrentImageNameWithoutExtension() 
String savePath = buildFilePath(saveDir, "${saveName}.zip")
exportAllObjectsToGeoJson(savePath)