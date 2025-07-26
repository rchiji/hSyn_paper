import qupath.lib.regions.ImagePlane
import ij.IJ
import ij.process.ColorProcessor
import qupath.imagej.processing.RoiLabeling
import qupath.imagej.tools.IJTools

// ----- parameter setting -------------------------------------------------

// マスク画像が入っているフォルダへのパス
def dirLabel = "C:/Users/admin/Documents/東京大学 整形外科/ヒト滑膜/hSyn_QuPath3/Train_add3"
def extension = ".png"

List className = ["Background","Immune cells", "plasma", "Fibro(loose)", "Fibro(dense,regular)",
"Fibro(dense,irregular)", "lining", "vessel", "vessel(large)",
"adipose","Stroma", "muscle", "RBC"]
//def className = null
// ファイル名にd=でダウンサンプリング係数が入っていない場合はここで設定。
def downsample = 1 // nullだとd=のところ探しに行く
// ---------------------------------------------------------------------------


// スライド名を取得
slidename = GeneralTools.getNameWithoutExtension(getCurrentImageData().getServer().getMetadata().getName())

// ローカルからスライド名をファイル名に含む画像一覧を取得
File folder = new File(dirLabel)
File[] listOfFiles = folder.listFiles({f -> f.getName().contains(slidename) && f.getName().endsWith(extension)} as FileFilter) as List
print listOfFiles

// 1画像ずつ処理 (もし複数枚あった時用にeach文を組んでいる。)
listOfFiles.each { file ->
    def path = file.getPath()
    
    // ファイル名のd=のところからダウンサンプリング係数を調べる。
    if (downsample == null){
        def filename = file.getName()
        def parts = filename.split(" ")
        def regionParts = parts[-1].split(",") as List;
        if (regionParts.size() == 4) {
            regionParts[0] = regionParts[0][1..-1]
            regionParts.add(0, "[d=1")
        }
        downsample = regionParts[0].replace("[d=", "") as float
        print "Downsampling factor is " + downsample
        }
    
    def imp = IJ.openImage(path) // ImageJで画像を読み込み

    print "Now processing: " + path
   
    ImagePlane plane = ImagePlane.getDefaultPlane()
    
    // Convert labels to ImageJ ROIs
    def ip = imp.getProcessor()
       // --> 読み込んだ画像のこんな情報 ip[width=6968, height=8048, bits=8, min=0.0, max=255.0]
    if (ip instanceof ColorProcessor) {
        throw new IllegalArgumentException("RGB images are not supported!")
    }
    // 最大輝度値から分類クラス数を出す。
    int n = imp.getStatistics().max as int
    print n + " Classes are identified"
    if (n == 0) {
        print 'No objects found!'
        return
    }
    
    // それぞれのclassのROIに変換
    def roisIJ = RoiLabeling.labelsToConnectedROIs(ip, n)
    print roisIJ
    
    // Convert ImageJ ROIs to QuPath ROIs
    def rois = roisIJ.collect {
        if (it == null)
            return
        return IJTools.convertToROI(it, 0, 0, downsample, plane);
    }
    
    for(i=0; i<=n; i++){   
        def roi = rois[i]
        print roi
        if (roi == null){
            continue // そのクラスのroiが無かったらfor文をskipして次の繰り返しに進む。
        }
            
        
        // Convert QuPath ROIs to objects
        def pathObject = PathObjects.createAnnotationObject(roi)
        if (className == null){
            addObjects(pathObject)
        } else {  
            pathObject.setPathClass(getPathClass(className[i]))
            addObjects(pathObject)
        }
        }
        }
