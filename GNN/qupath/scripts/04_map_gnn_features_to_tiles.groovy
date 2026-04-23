// Parameters ----------------------------------------------------
String path = buildFilePath(PROJECT_BASE_DIR, "../data/SLICTile_Features_241125",
    getCurrentImageName().split("_HE")[0] + ".txt")

String delim = "\t"
def objects = getTileObjects()
// ---------------------------------------------------------------

import qupath.lib.gui.measure.ObservableMeasurementTableData

def lines = new File(path).readLines()

def header = lines.pop().split(delim)

Map objectsById = objects.groupBy(object -> object.getID().toString())

for (def line in lines) {
    Map map = [:]
    def content = line.split(delim)
    for (int i = 0; i < header.size(); i++) {
        if ( header[i] != "ID") {
            map[header[i]] = content[i] as double
        } else {
            map[header[i]] = content[i]
        }
        }
    
    String id = map['ID']
    def object = objectsById[id][0]
    map.remove("ID")
    
    if (object == null) {
        continue
    }
    def ml = object.getMeasurementList()
    ml.putAll(map)
    }