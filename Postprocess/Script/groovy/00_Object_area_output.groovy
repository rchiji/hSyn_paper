import qupath.lib.gui.measure.ObservableMeasurementTableData
// import qupath.lib.gui.commands.SummaryMeasurementTableCommand

// Objects to export
List objects = getAnnotationObjects()
// Get the list of MeasurementList entries

// Create an instance of ObservableMeasurementTableData
def ob = new ObservableMeasurementTableData()
// Register objects to export
ob.setImageData(getCurrentImageData(), objects)
List terms = ob.getAllNames()

// Path to the output directory
String saveDir = buildFilePath(PROJECT_BASE_DIR, "measure/Annotation_ratio")
mkdirs saveDir
// Output file for measurement results
String imageName = getCurrentImageNameWithoutExtension()
String filePath = buildFilePath(saveDir, imageName + ".txt")

// Create a File object
File file = new File(filePath)

// Delete the file if it already exists
file.delete()  // Remove this line to append instead
// Create a new file
file.createNewFile()


// Use StringBuilder to store content in a buffer
StringBuilder buffer = new StringBuilder()

// Write header
buffer.append("Object ID" +"\t")
terms.eachWithIndex { key, i ->
    buffer.append(key)
    if ( i < terms.size() -1 ) {
        buffer.append("\t") 
        } else {
        buffer.append("\n")
        }
    }

objects.each {
    // Get Object ID to export
    String id = it.getID()
    buffer.append( "${id}\t" ) // append
    
    // Get measurement values
    def ml = it.getMeasurementList()
    
    // Extract and append measurement items
    terms.eachWithIndex { key, i ->
        def value = ob.getNumericValue(it, key)
        if ( value.isNaN() ) {
            value = ob.getStringValue(it, key)
            }
        
        buffer.append( value )
        if ( i < terms.size() -1 ) { 
            buffer.append("\t")
            } else {
            buffer.append("\n")
            }
        }
    }



file.withWriter { writer ->
    writer.write( buffer.toString() )
    }