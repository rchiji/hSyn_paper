//-------------------------
List classNames = ["Background","Immune cells", "plasma", "Fibro(loose)", "Fibro(dense,regular)",
"Fibro(dense,irregular)", "lining", "vessel", "vessel(large)",
"adipose", "nerve", "Stroma", "muscle", "RBC"]
def downsample = 5
def extension = ".tif"
def mask_extension = ".png"
def saveDir = "Train_add3"
//--------------------------

def imageData = getCurrentImageData()
def imageName = GeneralTools.getNameWithoutExtension(imageData.getServer().getMetadata().getName())
//Make sure the location you want to save the files to exists - requires a Project
def pathOutput = buildFilePath(PROJECT_BASE_DIR, saveDir)
mkdirs(pathOutput)

// Create an ImageServer where the pixels are derived from annotations
def labelServer = new LabeledImageServer.Builder(imageData)
  .backgroundLabel(0, ColorTools.WHITE) // Specify background label (usually 0 or 255)
  .downsample(downsample)
  .multichannelOutput(false)

def counter = 1
classNames.each {
    labelServer.addLabel(it, counter)
    counter++
    }
labelServer = labelServer.build()

def annotations = getSelectedObjects()
annotations.each{anno -> 
    def roi = anno.getROI()
    int h = roi.getBoundsHeight()
    int w = roi.getBoundsWidth()
    int x = roi.getBoundsX()
    int y = roi.getBoundsY()
    
    def tileInfo = "[d="+downsample + ",x="+x+",y="+y+",w="+w+",h="+h+"]"
    def saveName = imageName + " " + tileInfo + extension
    print saveName
    def savePath = buildFilePath(pathOutput, saveName)    
    def requestROI = RegionRequest.createInstance(getCurrentServer().getPath(), downsample, roi)
    writeImageRegion(getCurrentServer(), requestROI, savePath)
    
    // マスクの書き出し
    def mask_saveName = imageName + " " + tileInfo + mask_extension
    def mask_savePath = buildFilePath(pathOutput, mask_saveName)
    def requestMaskROI = RegionRequest.createInstance(labelServer.getPath(), downsample, roi) // label serverからROIをリクエスト    
    writeImageRegion(labelServer, requestROI, mask_savePath)
    
    removeObject(anno, true)
}


