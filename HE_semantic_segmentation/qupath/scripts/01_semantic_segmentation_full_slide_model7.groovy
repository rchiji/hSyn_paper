// ----- PARAMETERS --------------
// Downsampling factor
int downsample = 5
// Pixel size of each tile edge
int tileSize = 512
// Overlap size between tiles (in pixels)
int overlap = 256

def imageExtension = ".jpg"

def targets = null // getSelectedObject()

// List of PathClass for prediction targets (multiple classes can be specified) e.g., ["Region*"]
def targetClass = null // If null, the entire image is used
// Remove annotations other than the specified targetClass (if any)
boolean deleteAnnos = true
// ----------------------------------------

// Path to python.exe
String pythonEnvPath = "C:/Users/admin/anaconda3/envs/tf/python.exe"
// Path to the Python script
String pythonFile = "scripts/area_prediction.py"

// Path to the model
String modelPath = "../model_weights/model7_241111.h5"
// custom_scriptsパッケージのパス
String packagePath =  "scripts/custom_scripts"

// Batch size
int batch_size = 100
boolean centering = true
boolean average_probability = true
int label_splitsize = -1
// ----------------------------------------

List classNames = [
    "Immune cells", 
    "plasma",
    "Fibro(loose)",
    "Fibro(dense,regular)",
    "Fibro(dense,irregular)",
    "lining", 
    "vessel",
    "vessel(large)",
    "adipose",
    "Stroma",
    "muscle",
    "RBC"
    ]
// --------------------------------

String slideName = getCurrentImageNameWithoutExtension()
def pathOutput = new File(buildFilePath(PROJECT_BASE_DIR, "tmp", slideName))
String outPath = pathOutput.getAbsolutePath()
mkdirs(outPath)


def label_dir = new File(buildFilePath(PROJECT_BASE_DIR, "wsi_labels_model7_241111", slideName))
label_dir = label_dir.getAbsolutePath()


    
print "1. Tile exported"
tileExport(targets, targetClass, deleteAnnos, downsample,
           imageExtension, tileSize, overlap, outPath)
releaseMemory()    
    
print "2. Python process" 
def command = [
    pythonEnvPath,
    pythonFile,
    modelPath,
    packagePath,
    outPath,
    "--batch_size=${batch_size}",
    "--tile_size=${tileSize}",
    "--save_dir=${label_dir}",
    "--centering=${centering}",
    "--average_probability=${average_probability}",
    "--label_splitsize=${label_splitsize}"
    ]
    
pythonRunner(command)
releaseMemory()


print "3. Importing QuPath Objects"
importObjects(label_dir, slideName, downsample, classNames)
releaseMemory()

 
print "4. Delete tmp dir"
pathOutput.deleteDir()
clearSelectedObjects()
print "Done!!"



def releaseMemory() {
    Thread.sleep(100)
    javafx.application.Platform.runLater {
        getCurrentViewer().getImageRegionStore().cache.clear()
        System.gc()
    }
    Thread.sleep(100)    
    }
    
def tileExport(targets, targetClass, deleteAnnos, downsample,
               imageExtension, tileSize, overlap, outPath) {
    def imageData = getCurrentImageData()
    
    boolean annotatedTilesOnly = true
    if (targets == null) {
        if (targetClass == null) {
            annotatedTilesOnly = false
            } else {
            targets = getAnnotationObjects().findAll {
                it.getPathClass().toString() in targetClass
            }
            }
        }
    
    if (targets != null) {
        if(deleteAnnos == true) {
            del = getAnnotationObjects() - targets
            removeObjects(del, false)
            }    
        
        def allObjects = getAllObjects(false)
        removeList = allObjects - targets
        if( removeList.size() > 0) {
            removeObjects(removeList,true)
            }    
        } else {
            removeList = []    
        }
    
    new TileExporter(imageData)
        .downsample(downsample)
        .imageExtension(imageExtension)
        .tileSize(tileSize)
        .overlap(overlap)
        .annotatedTilesOnly(annotatedTilesOnly)
        .includePartialTiles(true)
        .writeTiles(outPath)
    
    addObjects(removeList)    
    }


def pythonRunner(command, boolean gpu=true) {   
    def pb = new ProcessBuilder(*command).redirectErrorStream( true )
        
    File tfPath = new File(command[0]).parentFile
    File sysPath = new File(tfPath, "Library\\bin")
    
    if (sysPath != null & gpu == true) {
        def env = pb.environment()
        def existingPath = env['Path']
        def updatedPath = existingPath + File.pathSeparator + sysPath.toString()
        env['Path'] = updatedPath         
        }        
        
    def out = new StringBuilder()    
    def process = pb.start()
    process.consumeProcessOutput(out, out)
    
    while ( process.isAlive() ) {
        if ( out.size() > 0 ) {
            logger.info( out.toString() )
            out.setLength( 0 )
        }
        sleep(200)
        }
    }

import ij.IJ
import qupath.imagej.processing.RoiLabeling

def importObjects(label_dir, slideName, downsample, classNames) {
    ImagePlane plane = ImagePlane.getDefaultPlane()
    
    File folder = new File(label_dir)
    File[] listOfFiles = folder.listFiles(
        {f -> f.getName().contains(slideName) && f.getName().endsWith("png")} as FileFilter
        ) as List
    
    Arrays.stream(listOfFiles).parallel().forEach { file ->
        def path = file.getAbsolutePath()  
        def parts = file.getName().split(" ")
        def regionParts = parts[-1].split(",") as List
        def xPosition = regionParts[1].replace("x=", "") as int
        def yPosition = regionParts[2].replace("y=","") as int
        def imp = IJ.openImage(path) 
        def ip = imp.getProcessor()
        int n = imp.getStatistics().max as int 
        if(n == 0){return}
        
        List roisIJ = RoiLabeling.labelsToConnectedROIs(ip, n)
        
        List pathObjects = []
    
        roisIJ.parallelStream().forEach { 
            if (it == null) {
                return
                }
            int idx = roisIJ.indexOf(it)    
            def ROI = IJTools.convertToROI(it,-xPosition/downsample,-yPosition/downsample,downsample,ImagePlane.getDefaultPlane())
            def pathObject = PathObjects.createAnnotationObject(ROI, getPathClass(classNames[idx]))
            pathObjects << pathObject
            }
        
        pathObjects = pathObjects.findAll{it != null}
        
        if(pathObjects.size() > 0) {
            addObjects(pathObjects)
        }
    }    
    }