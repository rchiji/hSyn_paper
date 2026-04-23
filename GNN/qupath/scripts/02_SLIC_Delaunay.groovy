removeObjects( getTileObjects(), true)

List objects = getAnnotationObjects()
def mergeObject = PathObjectTools.mergeObjects(objects)
mergeObject.setPathClass(getPathClass("Region*"))
addObjects(mergeObject)



selectObjectsByClassification('Region*')
runPlugin('qupath.imagej.superpixels.SLICSuperpixelsPlugin',\
    '{"sigmaMicrons":0.0,\
    "spacingMicrons":20.0,\
    "maxIterations":10,\
    "regularization":0.01,"\
    adaptRegularization":true}'
    )

