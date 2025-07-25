// Parameter ---------------
double downsample = 10
extension = ".jpg" // 保存時の拡張子
saveDir = "rendered_thumnails_241111" // 保存先のフォルダ
// -------------------------

// 追加で必要なmoduleをimport
import qupath.lib.gui.images.servers.RenderedImageServer
import qupath.lib.gui.viewer.overlays.HierarchyOverlay

// 保存時に使用するために画像名を取得
def server = getCurrentServer()
String imageName = getCurrentImageNameWithoutExtension()

// 保存先のフォルダを作成
String pathOutput = buildFilePath(PROJECT_BASE_DIR, saveDir)
mkdirs(pathOutput)

// objectの表示設定を取得
def viewer = getCurrentViewer()
def options = viewer.getOverlayOptions()

// RenderedImageServerを構築
def imageData = getCurrentImageData()
def renderedServer = new RenderedImageServer.Builder(imageData)
    .downsamples(downsample)
    .layers(new HierarchyOverlay(null, options, imageData))
    .build()

// renderedServerから書き出し
String savePath = "${pathOutput}/${imageName}${extension}"
writeImage(renderedServer, savePath)