List<PathObject> objects = getAnnotationObjects()

def regionObject = PathObjectTools.mergeObjects( objects )
regionObject.setPathClass( getPathClass("Region*") )
addObject( regionObject )