// ----- Tile作成PARAMETERS --------------
// ダウンサンプリング係数
int downsample = 5
// タイル辺のピクセルサイズ
int tileSize = 512
// 重複タイルピクセル数
int overlap = 256
// 拡張子
def imageExtension = ".jpg"

// 予測対象のObject
def targets = null // getSelectedObject()

// 予測対象のPathClassリスト（複数class指定可能） Ex) ["Region*"]
def targetClass = null // nullなら画像全体
// 既にあるtargetClass以外のアノテーションを削除
boolean deleteAnnos = true
// ----------------------------------------

// ----- Python実行PARAMETERS --------------
// python.exeまでのパス
String pythonEnvPath = "C:\\Users\\admin\\anaconda3\\envs\\tf\\python.exe"
// pythonスクリプトへのパス
String pythonFile = "C:\\Users\\admin\\Documents\\東京大学 整形外科\\ヒト滑膜\\hSyn_QuPath3\\scripts\\area_prediction.py"

// modelへのパス
String modelPath = "C:\\Users\\admin\\Documents\\東京大学 整形外科\\ヒト滑膜\\hSyn_QuPath3\\model7_241111.h5"
// custom_scriptsパッケージのパス
String packagePath =  "C:/Users/admin/Documents/東京大学 整形外科/ヒト滑膜/hSyn_QuPath3/scripts/custom_scripts"

// バッチサイズ
int batch_size = 100
boolean centering = true
boolean average_probability = true
int label_splitsize = -1
// ----------------------------------------

// ----- Object importのPARAMETERS ---------
// アノテーションクラス名。輝度値順と揃えておく必要アリ。
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
// 出力先のフォルダ作成
def pathOutput = new File(buildFilePath(PROJECT_BASE_DIR, "tmp", slideName))
String outPath = pathOutput.getAbsolutePath()
mkdirs(outPath)


// 予測後ラベル画像の保存場所
def label_dir = new File(buildFilePath(PROJECT_BASE_DIR, "wsi_labels_model7_241111", slideName))
label_dir = label_dir.getAbsolutePath()


    
// 1. 予測対象領域のタイル書き出し
print "1. Tile exported"
tileExport(targets, targetClass, deleteAnnos, downsample,
           imageExtension, tileSize, overlap, outPath)
releaseMemory()    
    
// 2. Pythonスクリプトの実行     
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

// 3. QuPath Object作成  

print "3. Importing QuPath Objects"
importObjects(label_dir, slideName, downsample, classNames)
releaseMemory()

// 4. tmp folderの削除
 
print "4. Delete tmp dir"
pathOutput.deleteDir()
clearSelectedObjects()
print "Done!!"



def releaseMemory() {
    // CPUメモリ解放
    Thread.sleep(100)
    // Try to reclaim whatever memory we can, including emptying the tile cache
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
    // targetsが事前に無ければ、targetClassのObjectが予測対象となる
    if (targets == null) {
        if (targetClass == null) {
            // 対象クラスの指定も無ければ画像全体が対象
            annotatedTilesOnly = false
            } else {
            targets = getAnnotationObjects().findAll {
                it.getPathClass().toString() in targetClass
            }
            }
        }
    
    // targets以外のObjectは一旦削除
    if (targets != null) {
        // 既にあるtargets以外のアノテーションを削除
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
    
    // タイル画像の書き出し
    new TileExporter(imageData)
        .downsample(downsample)          // ダウンサンプリング係数
        .imageExtension(imageExtension)  // 画像の書き出し拡張子
        .tileSize(tileSize)             // タイル辺のピクセル数
        .overlap(overlap)               // 隣り合うタイルの重複幅
        .annotatedTilesOnly(annotatedTilesOnly)   // trueにするとタイル画像中にラベルがあるタイルのみが書き出される
        .includePartialTiles(true)  // タイルサイズに満たない画像を書き出すかどうか
        .writeTiles(outPath)  // タイル画像の書き出し先のパス
    
    // Objectを元に戻す
    addObjects(removeList)    
    }


def pythonRunner(command, boolean gpu=true) {   
    // プロセスを作成
    def pb = new ProcessBuilder(*command).redirectErrorStream( true )
        
    // Fileクラスを使用してパスを構築
    File tfPath = new File(command[0]).parentFile
    File sysPath = new File(tfPath, "Library\\bin")
    
    // ---- システム環境変数にcudaがあるbinフォルダのパスを追記 ----
    if (sysPath != null & gpu == true) {
        def env = pb.environment()
        // 既存のPATHを取得
        def existingPath = env['Path']
        // 既存のPATHに新しいパスを追加
        def updatedPath = existingPath + File.pathSeparator + sysPath.toString()
        // プロセス内でPATHを更新
        env['Path'] = updatedPath         
        }        
        
    def out = new StringBuilder()    
    // プロセスを実行
    def process = pb.start()
    // プロセスの出力を受け取る
    process.consumeProcessOutput(out, out)
    
    // プロセス実行中
    while ( process.isAlive() ) {
        if ( out.size() > 0 ) {
            logger.info( out.toString() ) // print out.toString()でも同じ
            out.setLength( 0 )
        }
        sleep(200)
        }
    }

import ij.IJ
import qupath.imagej.processing.RoiLabeling

def importObjects(label_dir, slideName, downsample, classNames) {
    ImagePlane plane = ImagePlane.getDefaultPlane()
    
    // スライド名をファイル名に含むラベル画像一覧を取得
    File folder = new File(label_dir)
    File[] listOfFiles = folder.listFiles(
        {f -> f.getName().contains(slideName) && f.getName().endsWith("png")} as FileFilter
        ) as List
    
    // 画像中の輝度値ごとにImageJ ROIを作成 -> QuPathアノテーションに変換
    Arrays.stream(listOfFiles).parallel().forEach { file ->
        def path = file.getAbsolutePath() // マスク画像のpath取得  
        // ファイル名からタイルの貼り付け位置を取得。
        def parts = file.getName().split(" ")
        def regionParts = parts[-1].split(",") as List
        // def downsample = regionParts[0].replace("[d=", "") as int
        def xPosition = regionParts[1].replace("x=", "") as int
        def yPosition = regionParts[2].replace("y=","") as int
        // ImageJでマスク画像を開く
        def imp = IJ.openImage(path) 
        // ImageProcessorを取得
        def ip = imp.getProcessor()
        // 最大輝度値から何クラスあるか記録
        int n = imp.getStatistics().max as int 
        //print n + " Classes are identified"
        if(n == 0){return}
        
        // ラベル画像からImageJ Roiを作成
        List roisIJ = RoiLabeling.labelsToConnectedROIs(ip, n)
        
        // ImageJ RoiをQuPath Objectに変換
        List pathObjects = []
    
        roisIJ.parallelStream().forEach { 
            // 輝度値iのROIが無ければreturn。返り値はnull
            if (it == null) {
                return
                }
            int idx = roisIJ.indexOf(it)    
            // 輝度値iのROIがあればQuPath ROIに変換し、さらにQuPath Objectに変換。i番目のPathClassをつける。
            def ROI = IJTools.convertToROI(it,-xPosition/downsample,-yPosition/downsample,downsample,ImagePlane.getDefaultPlane())
            def pathObject = PathObjects.createAnnotationObject(ROI, getPathClass(classNames[idx]))
            pathObjects << pathObject
            }
        
        // リストからnullを除く
        pathObjects = pathObjects.findAll{it != null}
        
        // QuPathに反映
        if(pathObjects.size() > 0) {
            addObjects(pathObjects)
        }
    }    
    }