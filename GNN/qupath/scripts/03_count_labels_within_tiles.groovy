// Parameters ------
double downsample = 1
def objects = getTileObjects()

String saveDir = "../data/labelRatio_241117"
mkdirs buildFilePath(PROJECT_BASE_DIR, saveDir)

String imageName = getCurrentImageNameWithoutExtension()
String filePath = buildFilePath(PROJECT_BASE_DIR, saveDir, imageName + ".txt")

List className = ["Background","Immune cells", "plasma", "Fibro(loose)", "Fibro(dense,regular)",
"Fibro(dense,irregular)", "lining", "vessel", "vessel(large)",
"adipose", "Stroma", "muscle", "RBC"]
// -----------------

import ij.ImagePlus
import ij.IJ

def server = getCurrentServer()

List labelTables = []
objects.each {
    def roi = it.getROI()
    
    int xPosition = roi.getBoundsX()
    int yPosition = roi.getBoundsY()
    
    def request = RegionRequest.createInstance(server.getPath(), downsample, roi)
    
    ImagePlus imp = IJTools.convertToImagePlus(server, request).getImage()
    
    int width = imp.getWidth()
    int height = imp.getHeight()

    def roiIJ = IJTools.convertToIJRoi(roi,-roi.getBoundsX(),-roi.getBoundsY(),1)
    
    imp.setRoi(roiIJ)
    
    Map labelTable = [:]
    labelTable["ID"] = it.getID()
    labelTable["x"] = roi.getCentroidX()
    labelTable["y"] = roi.getCentroidY()
    
    (1..12).each { label ->
        labelTable.put(label,0)
        }
    
    for( int x in 0..width-1 ) {
       for( int y in 0..height-1 ) {
           if ( roiIJ.contains(x,y) ) {
               label = imp.getPixel(x,y)
               if (label[0] > 0) {
                   labelTable[label[0]] += 1               
                   } // <- if (label[0] > 0)               
               } // <- if ( roiIJ.contains(x,y) )
           } // <- for( int y in 0..height-1 )
       } // <- for( int x in 0..width-1 )
    
    labelTables << labelTable
    
    List<Number> values = labelTable.values().toArray()[3..(labelTable.size()-1)]
    int maxValue = values.max()
    int key = labelTable.find{ key,value ->
        value == maxValue
        }.key

    it.setName(className[key])
    def ml = it.getMeasurementList()
    ml.put("Max label: ",key)
    
    values.eachWithIndex{ v, i ->
        ml.put("Px num: " + className[i+1], v)
        }
    }
    
    
print labelTables

print "Writing..."

new File(filePath).withWriter { writer ->
    def bufferedWriter = new BufferedWriter(writer)
    keys = labelTables[0].keySet()
    for(String key:keys) {
        bufferedWriter.write(key + "\t")
        }
    bufferedWriter.write("\n")
    
    StringBuilder sb = new StringBuilder()
    labelTables.each { labelTable ->
        Set<String> keys = labelTable.keySet()
        
        sb.setLength(0)  // Clear the StringBuilder
        for(key:keys) {
            sb.append(labelTable.get(key))
            sb.append("\t")
            }
        sb.append('\n')
        bufferedWriter.write(sb.toString())
    }
    bufferedWriter.flush()
}

print "Done!!"    