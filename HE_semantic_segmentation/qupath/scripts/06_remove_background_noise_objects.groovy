// --- 1. Split annotation objects ---
selectAnnotations();
runPlugin('qupath.lib.plugins.objects.SplitAnnotationsPlugin','{}')
resetSelection()

// --- 2. Create Ignore* objects in background regions ---

// Parameter -------------------------------------------
def op = ImageOps.buildImageDataOp()
                 .appendOps(
                     ImageOps.Channels.mean(),
                     ImageOps.Filters.gaussianBlur(5),
                     ImageOps.Threshold.threshold(220)
                     )
                

double downsample = 64

Map classifications = Map.of(
    1, getPathClass("Ignore*")
    )

double minArea = 10000000
double minHoleArea = 0
boolean INCLUDE_IGNORED = true
boolean DELETE_EXISTING = false
boolean SELECT_NEW = false
boolean SPLIT = false
// -----------------------------------------------------

import qupath.opencv.ml.pixel.PixelClassifiers

def imageData = getCurrentImageData()

def resolution = imageData.getServer().getPixelCalibration().createScaledInstance(downsample,downsample)

def classifier = PixelClassifiers.createClassifier(op, resolution, classifications)


List options = [
    INCLUDE_IGNORED == true ? PixelClassifierTools.CreateObjectOptions.INCLUDE_IGNORED : null,
    DELETE_EXISTING == true ? PixelClassifierTools.CreateObjectOptions.DELETE_EXISTING : null,
    SELECT_NEW == true ? PixelClassifierTools.CreateObjectOptions.SELECT_NEW : null,
    SPLIT == true ? PixelClassifierTools.CreateObjectOptions.SPLIT : null    
    ]

PixelClassifierTools.createAnnotationsFromPixelClassifier(
    imageData,
    classifier,
    minArea,
    minHoleArea,
    *options)


// --- 3. Remove objects contained within Ignore* ---
def ignoreROI = getAnnotationObjects().find {
    it.getPathClass().toString() == "Ignore*"
    }.getROI()

def hierarchy = getCurrentHierarchy()
List deleteObjects = hierarchy.getObjectsForROI( null, ignoreROI )

removeObjects(deleteObjects, true)

def ignores = getAnnotationObjects().findAll {
    it.getPathClass().toString() == "Ignore*"
    }
if( ignores.size() > 0) {
    removeObjects( ignores, true)
    }
