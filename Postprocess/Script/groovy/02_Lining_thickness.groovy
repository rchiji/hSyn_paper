// ----- Parameters -----------------
// sample name
def SampleName = getProjectEntry()?.getImageName()?.replaceFirst(/\.[^.]+$/, '') ?: "Unknown"

// Path
def outDir = buildFilePath(PROJECT_BASE_DIR, "measure/thickness_lining")
def outFile = new File(outDir, "${SampleName}_thickness_lining.txt")

outFile.getParentFile().mkdirs()

// Object
def targetClass = getPathClass('lining')
def targetObjects = getAnnotationObjects().findAll {
    it.getPathClass() == targetClass
}
if (targetObjects.isEmpty()) {
    print "No 'lining' objects found in this image."
    return 
}

// Parameter
double downsample = 3
boolean prune = false
def minDiameter = 20
boolean calibration = true
boolean showSkelton = true
boolean showSkeltonPoint = true
boolean showCircle = true
// -----------------------------------

// remove Circle, Skelton class
targetObjects = targetObjects - targetObjects.findAll {
    it.getPathClass().toString() in ["Circle","Skelton"]
    } 

if(calibration) {
    pixelSize = getCurrentServer().getMetadata().getAveragedPixelSize()
    } else {
    pixelSize = 1
    }
print "Pixel size: ${pixelSize}/px"

targetObjects.each { targetObj ->
    removeObjects(targetObj.getChildObjects().findAll {
        it.getPathClass().toString() in ["Circle","Skelton"]
        },false)

    def targetROI = targetObj.getROI()
    
    if( downsample < 0 ) {
        downsample = estimateDownsample(targetROI)
        print "Downsample: ${downsample}"
    }    
    
    def bbROI = getBBROI(targetROI)
    def outerROI = RoiTools.subtract(bbROI, targetROI)
    
    List skeletonCoords = getSkeletonCoord(targetROI,downsample,prune,showSkelton)
    
    if( showSkeltonPoint ) {
        addObjects(PathObjects.createAnnotationObject(
            ROIs.createPointsROI(skeletonCoords, null), getPathClass("Skelton"))
            )        
        }
        
    List results = skeletonCoords.collect {
        def pointROI = ROIs.createPointsROI([it], null)
        double diameter = RoiTools.getBoundaryDistance(pointROI, outerROI) * 2
        return [it, diameter]
        }
    
    // ----- Statistics -----
    if (minDiameter != null) {
        results = results.findAll { it[1] > minDiameter }
    }

    if (results.isEmpty()) {
        print "No valid diameter points after filtering; skip measurements for this object."
        return
    }

    List<Double> diameters = results.collect { it[1] * pixelSize as double }

    double mean = diameters.average()
    double min  = diameters.min()
    double max  = diameters.max()

    double sd
    if (diameters.size() < 2) {
        sd = 0.0
        print "Only 1 diameter point; set SD=0."
    } else {
        sd = calcSD(diameters)
    }

    def ml = targetObj.getMeasurementList()
    ml.put("Diameter mean", mean)
    ml.put("Diameter min",  min)
    ml.put("Diameter max",  max)
    ml.put("Diameter SD",   sd)
    ml.close()

    print "Mean: ${mean}"
    print "Min: ${min}"
    print "Max: ${max}"
    print "SD: ${sd}"
    
    if (showCircle) {
        if (results.isEmpty()) {
            print "No skeleton points; skip circles."
        } else {
            Random rnd = new Random()
            // 0 から results.size()-1 まで
            List tmp = (0..<results.size()).toList()
            int nDraw = Math.min(100, tmp.size()) 
            def indexes = []
            for (int k = 0; k < nDraw; k++) {
                int j = rnd.nextInt(tmp.size())
                indexes.add(tmp.remove(j))
            }
            List circleROIs = []
            for (index in indexes) {
                def res = results[index]
                def circleROI = createCircleROI(
                    res[0].getX(),
                    res[0].getY(),
                    res[1]
                )
                circleROIs << circleROI
            }
            def circleObjs = circleROIs.collect {
                PathObjects.createDetectionObject(
                    it,
                    getPathClass("Circle")
                )
            }
            addObjects(circleObjs)
        }
    }
}


// Save measurements
try {
    this.&saveAnnotationMeasurements(outFile)
    print "Saved (QuPath export): " + outFile
} catch (Exception e) {

    outFile.withWriter('UTF-8') { w ->
        w.writeLine("Sample\tObjectID\tClass\tDiameter mean\tDiameter min\tDiameter max\tDiameter SD")

        targetObjects.each { obj ->
            def ml = obj.getMeasurementList()
            def mean = ml.get("Diameter mean")
            def min  = ml.get("Diameter min")
            def max  = ml.get("Diameter max")
            def sd   = ml.get("Diameter SD")

            w.writeLine([
                SampleName,
                obj.getID(),
                obj.getPathClass()?.toString(),
                mean, min, max, sd
            ].join('\t'))
        }
    }
    print "Saved (TSV fallback): " + outFile
}



/*
 * Bounding box
 */
import qupath.lib.geom.Point2
def getBBROI(roi,dilate=1) {
    def x = roi.getBoundsX() - dilate
    def y = roi.getBoundsY() - dilate
    def w = roi.getBoundsWidth() + dilate*2
    def h = roi.getBoundsHeight() + dilate*2
    
    List bbCoord = [
        new Point2(x,y), 
        new Point2(x+w,y),
        new Point2(x+w,y+h),
        new Point2(x,y+h)
        ]
        
    def BBROI = ROIs.createPolygonROI(bbCoord,null)
    
    return BBROI
    }


/*
 * ImageJ Skeleton ➔ Extract skeletal coordinates.
 */
import ij.ImagePlus
import ij.IJ

def getSkeletonCoord(roi,downsample,prune,addObject) {
    
    def imgMask = BufferedImageTools.createROIMask(roi, downsample)
    def mask = new ImagePlus("Mask", imgMask)
    
    IJ.run(mask, "Skeletonize (2D/3D)", "");
    if( prune) {
        IJ.run(mask, "Analyze Skeleton (2D/3D)", "prune=none prune_0");
    }
    
    IJ.setRawThreshold(mask, 1, 255); 
    IJ.run(mask, "Convert to Mask", ""); 
    
    int width = mask.getWidth()
    int height = mask.getHeight()
    List truePixels = []
    for( int x in 0..width ) {
       for( int y in 0..height ) {
           if( mask.getPixel(x,y)[0] > 0 ) {
               x2 = x * downsample + roi.getBoundsX() 
               y2 = y * downsample + roi.getBoundsY() 
               truePixels << new Point2(x2,y2)
           }
       }
    }
    
    if( addObject) {
        IJ.run(mask, "Create Selection", ""); 
        def newRoiIJ = mask.getRoi()
        def newROI = IJTools.convertToROI(newRoiIJ, -roi.getBoundsX()/downsample, -roi.getBoundsY()/downsample, downsample, null)
        def pathObject = PathObjects.createAnnotationObject(newROI, getPathClass("Skelton"))
        addObjects(pathObject)
        }
        
    return truePixels
    }
    
/*
 * A function that generates a perfect circle given the center coordinates and diameter as arguments.
 */
def createCircleROI(centerX,centerY,diameter) {
    def roi = ROIs.createEllipseROI(0,0,diameter,diameter,null)
    def newROI = roi.translate(centerX - diameter/2, centerY - diameter/2)
    
    return newROI
    }    

/*
 * A function that estimates an appropriate downsampling factor based on the ROI area.
 */
def estimateDownsample(roi) {
    double area = roi.getArea()
    
    int digit = Math.floor(Math.log10(Math.abs(area)) + 1)
    
    double downsample = digit - 5
    if( downsample < 1) {
        downsample = 1 
    } else {
        downsample = downsample ** 2
    }
    return downsample
}

/*
 * A function for calculating the standard deviation.
 */    
def calcSD(list) {
    int n = list.size()
    if (n < 2) return 0.0   // or Double.NaN

    double mean = list.average()
    double variance = list.collect{ (it-mean)**2 }.sum() / (n - 1)
    double sd = Math.sqrt(variance)
    return Math.round(sd * 1e4) / 1e4
}

