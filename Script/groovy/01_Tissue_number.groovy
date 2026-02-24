// ----- Parameters -------------

// Class (The same applies to "Immune cells")
def targetClass = getPathClass("vessel(large)")

// Threshold ("Immune cells" -> 7000.)
double areaThreshold = 1500.0

def AreaKey = "Area µm^2"

// sample name
def SampleName = getProjectEntry()?.getImageName()?.replaceFirst(/\.[^.]+$/, '') ?: "Unknown"

// Path
def outDir = buildFilePath(PROJECT_BASE_DIR, "measure/count_vessel_large_filtered")
def outFile = new File(outDir, "${SampleName}_n_vessel_large.txt")

outFile.getParentFile().mkdirs()
// -------------------------------

// Annotation
def target_Annotations = getAnnotationObjects().findAll { it.getPathClass() == targetClass }

// Scale
def cal = getCurrentServer().getPixelCalibration()
double pixelArea = cal.getPixelWidthMicrons() * cal.getPixelHeightMicrons()

// Filtering
def filtered = target_Annotations.findAll { obj ->
    def roi = obj.getROI()
    if (roi == null) return false
    double area_um2 = roi.getArea() * pixelArea
    return area_um2 >= areaThreshold
}

// Count
def Count = filtered.size()

// Save
outFile.withWriter('UTF-8') { writer ->
    writer.writeLine("SampleName\tCount")
    writer.writeLine("${SampleName}\t${Count}")
}

print "Saved count for ${SampleName} to ${outFile.getAbsolutePath()}"
