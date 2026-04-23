// Parameter ---------------
double downsample = 10
extension = ".jpg"
saveDir = "predictions/rendered_thumnails_model7_241113"
// -------------------------

import qupath.lib.gui.images.servers.RenderedImageServer
import qupath.lib.gui.viewer.overlays.HierarchyOverlay

def server = getCurrentServer()
String imageName = getCurrentImageNameWithoutExtension()

String pathOutput = buildFilePath(PROJECT_BASE_DIR, saveDir)
mkdirs(pathOutput)

def viewer = getCurrentViewer()
def options = viewer.getOverlayOptions()

def imageData = getCurrentImageData()
def renderedServer = new RenderedImageServer.Builder(imageData)
    .downsamples(downsample)
    .layers(new HierarchyOverlay(null, options, imageData))
    .build()

String savePath = "${pathOutput}/${imageName}${extension}"
writeImage(renderedServer, savePath)