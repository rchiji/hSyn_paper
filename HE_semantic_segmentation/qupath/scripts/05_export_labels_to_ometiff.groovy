// Parameter ---------------
double downsample = 5
String extension = ".ome.tiff"
String saveDir = "predictions/label_ometiff_model7_241113"
// -------------------------

import qupath.lib.images.servers.LabeledImageServer

def imageData = getCurrentImageData()

String imageName = getCurrentImageNameWithoutExtension()
String pathOutput = buildFilePath(PROJECT_BASE_DIR, saveDir)
mkdirs(pathOutput)

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

def labelServer = new LabeledImageServer.Builder(imageData)
for(i = 0; i < classNames.size(); i++) {
    labelServer.addLabel(classNames[i], i+1)
    }
labelServer = labelServer
    .downsample(downsample)
    .multichannelOutput(false)
    .build()

def request = RegionRequest.createInstance(labelServer.getPath(), downsample,
                                            0,0, labelServer.getWidth(), labelServer.getHeight(), 0,0)

String savePath = "${pathOutput}/${imageName}${extension}"
writeImageRegion(labelServer, request, savePath) 

