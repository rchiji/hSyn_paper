String saveDir = buildFilePath(PROJECT_BASE_DIR, 'predictions/Prediction_label_json_model7_241113')
mkdirs(saveDir)

String saveName = getCurrentImageNameWithoutExtension() 
String savePath = buildFilePath(saveDir, "${saveName}.zip")
exportAllObjectsToGeoJson(savePath)