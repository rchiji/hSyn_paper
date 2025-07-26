// QuPathプロジェクトを取得
def project = getProject()

// QuPathプロジェクトの登録画像一覧を取得
List imgList = project.getImageList()

imgList.each{
    // 画像読み込み
    def img = it.readImageData()
    // アノテーション一覧を取得し、アノテーションが何もなければその画像はQuPathプロジェクトから削除
    if(img.getHierarchy().getAnnotationObjects().size==0){
        project.removeImage(it, true)
    }
}
