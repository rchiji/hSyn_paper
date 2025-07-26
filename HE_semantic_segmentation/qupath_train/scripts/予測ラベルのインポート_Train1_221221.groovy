// ----- SET THESE PARAMETERS -----
// マスク画像が入っているフォルダへのパス
def labelDir = buildFilePath(PROJECT_BASE_DIR, 'Prediction_label')
def extension = ".png"
// アノテーションクラス名。輝度値順と揃えておく必要アリ。
def className = ["Immune cells", "adipose", "vessel", "lining", "RBC", "Fibro(loose)",
"Fibro(Perivascular)","Fibro(dence,regular)","Fibro(dense,irregular)",
"cell infiltration","vessel(large)"]

// ダウンサンプリングしていれば要変更
def downsample = 5
// RoiManager
// https://javadoc.scijava.org/ImageJ1/ij/ij/plugin/frame/RoiManager.html
def saveDir = buildFilePath(PROJECT_BASE_DIR, 'Prediction_label_json')
// --------------------------------

import ij.plugin.frame.RoiManager
import ij.IJ
import qupath.imagej.tools.IJTools
ImagePlane plane = ImagePlane.getDefaultPlane()

// 保存先フォルダの作成
mkdirs(saveDir)

// スライド名を取得
slidename = GeneralTools.getNameWithoutExtension(getCurrentImageData().getServer().getMetadata().getName())

// ローカルからスライド名をファイル名に含む画像一覧を取得
File folder = new File(labelDir)
File[] listOfFiles = folder.listFiles({f -> f.getName().contains(slidename) && f.getName().endsWith(extension)} as FileFilter) as List

// 画像中の輝度値ごとにImageJ ROIを作成 -> QuPathアノテーションに変換
listOfFiles.each { file ->
    def path = file.getAbsolutePath() // マスク画像のpath取得  
    
    // ファイル名からタイルの貼り付け位置を取得。
    def parts = file.getName().split(" ")
    def regionParts = parts[-1].split(",") as List;
    def xPosition = regionParts[1].replace("x=", "") as int
    def yPosition = regionParts[2].replace("y=","") as int
    
    def imp = IJ.openImage(path) // ImageJでマスク画像を開く
    int n = imp.getStatistics().max as int // 最大輝度値から何クラスあるか記録
    print n + " Classes are identified"
    if(n == 0){return}
    
    def rm = new RoiManager() // ImageJ RoiManagerを起動
    
    // 輝度値ごとにループ。(輝度値0は背景なので、輝度値はi=1から開始)
    for(i=1; i<=n; i++){
        tmp = imp.duplicate(); // マスク画像を複製
        IJ.setThreshold(tmp, i, i); // 輝度値iでthresholding
        IJ.run(tmp, "Convert to Mask", ""); // 二値化
        IJ.run(tmp, "Create Selection", ""); // 選択
        
        // ROIを取得して、RoiManagerに登録　https://forum.image.sc/t/roimanager-cant-add-roi-if-identical-to-n-1/19298
        // ROIが無ければ次のfor loopへ
        if(tmp.getRoi()==null) {
            tmp.close(); continue 
        } else {
            rm.addRoi(tmp.getRoi())
        }
        
        rm.select(rm.getCount()-1); // RoiManagerの最後尾を選択。
        // QuPathアノテーションクラス名一覧を用意していなければ、index_輝度値の名前を付ける。
        // 用意していればclassNameリストのi-1番目の名前を付ける。
        if(className==null){
        rm.runCommand("Rename", "index_${i}");
        } else{
        rm.runCommand("Rename", className[i-1]);
        }
        tmp.close() // 複製画像は閉じる。
    }
    imp.close(); // マスク画像を閉じる。
    
    // RoiManagerのROIをQuPathアノテーションに変換
    for(j=0; j<rm.getCount();j++){
        // RoiManagerのj番目をROIクラスに変換
        def roi = IJTools.convertToROI(rm[j], -xPosition, -yPosition, downsample, plane)
        // ROIクラスをPathObjectsクラスに変換
        def pathObject = PathObjects.createAnnotationObject(roi)
        
        // ImageJ ROIのNameをアノテーションクラスに設定
        def pathClass = rm[j].getName()
        if (pathClass != null){
            pathObject.setPathClass(getPathClass(pathClass))
            }
        addObjects(pathObject) // QuPath上に反映
    }
    
    def saveName = GeneralTools.getNameWithoutExtension(file.getName()).split(" ")[0]
    def savePath = buildFilePath(saveDir, "${saveName}.zip")
    rm.runCommand("Save", savePath) // RoiManagerを保存
    //rm.reset() // RoiManagerの情報を消去
    rm.close() // RoiManagerを閉じる。
}

