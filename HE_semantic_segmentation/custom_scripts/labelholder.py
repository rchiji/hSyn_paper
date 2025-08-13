import os
import numpy as np
from tqdm.auto import tqdm
import scipy.ndimage
import cv2

from datasets import *
from tiling import *
from prediction import *

class LabelHolder:
    def __init__(self, path):
        self.path = path
        self.tile_info = getTileCoord(self.path)
        self.tile_info_d = downsampleCoord(self.tile_info)
        self.label = None
        self.label_final = None
        self.call_num = 0 # 何回呼ばれたか記録
        self.total_call_num = 0 # 合計何回呼ばれるかを事前計算した値
        self.surround = None
        self.region = None
        self.trimmed = False # crop済みかの記録
    
    def checkSurround(self):
        """
        周辺画像の予測ラベルが揃っているか確認
        """
        check = [surround.label is not None for surround in self.surround.values()]
        if check:
            return all(check)   
        else:
            return False
    
    def cropLabel(self,trim=None, overlap=None):
        """
        ラベル画像の上下左右からtrim幅分をトリミング
        """
        if trim is None:
            trim = overlap // 32
        if self.trimmed is False:
            self.label = self.label[trim:self.label.shape[0]-trim, trim:self.label.shape[1]-trim]
            self.trimmed = True
    
    def centerCrop(self, overlap):
        """
        ラベル画像の上下左右から重複幅/2をトリミング
        """
        trim = overlap // 2
        self.label_final = self.label_final[trim:self.label_final.shape[0]-trim, trim:self.label_final.shape[1]-trim]
        

    def averageProb(self, start_dict, num_classes, tile):
        """
        中心画像と周辺画像のoverlap領域のラベル確率を平均
        """
        tile_size = self.tile_info_d["w"]
        
        # 中心画像を登録
        y_start, x_start = start_dict[4]
        tile[y_start:y_start+self.label.shape[0], x_start:x_start+self.label.shape[1],:,4] = self.label
        
        # 周囲画像を登録
        for counter, slh in self.surround.items():
            y_start, x_start = start_dict[counter] 
            # labelを代入
            label = slh.label
            tile[y_start:y_start+label.shape[0], x_start:x_start+label.shape[1],:,counter] = label
            slh.call_num +=1
        
        # 中心画像の領域をcrop
        center = start_dict["center"]
        tile = tile[center:center+tile_size, center:center+tile_size]
        
        # 周囲画像の画像とで確率を平均 (tile_size,tile_size,num_classes,9) -> (tile,tile,num_classes)
        tile = np.einsum('ijkl->ijk', tile)/num_classes
        # 最大確率のラベルを採用 (tile_size,tile_size,num_classes) -> (tile_size,tile_size)
        tile = np.argmax(tile,axis=-1)
    
        self.label_final = tile
        
    def removeLabel(self):
        """
        画像がすでに必要回数呼び出されていたらラベル画像を削除
        """
        if self.call_num == self.total_call_num:
            del self.label


class LabelHolders:
    def __init__(self,
                 labelholders=None,                  
                 image_paths=None,
                 image_dir=None,
                ):
        
        # labelholderのリストを登録
        self.labelholders = labelholders
        self.image_paths = image_paths
        self.image_dir = image_dir
        if labelholders is None:           
            self.prepare_labelholders()
        
        self.filename = os.path.basename(self.labelholders[0].path)
        self.slidename = re.split("[ ]", self.filename)[0] # 半角スペースで区切ってスライド名を保存
        self.downsample = self.labelholders[0].tile_info["d"]
        self.tile_size = self.labelholders[0].tile_info_d["w"]
        self.num_classes = None
        # 隣接タイルとのoverlapや貼り付け座標のずれを計算
        self.calc_overlap()
        # whole slideサイズのゼロアレイ作成
        self.slide = createEmptySlide(labelholders=self.labelholders)
        # Center cropした場合の貼り付け先座標の更新
        self.register_croppedCoord()

    
    def prepare_labelholders(self):
        """
        明視野タイル画像のファイルパスリスト or ディレクトリパスから
        labelholderのリストを作成する
        """
        if self.image_paths is None:
            image_paths = prepare_image_paths(self.image_dir)
            
        image_paths = sorted(image_paths)
        self.labelholders = [ LabelHolder(path) for path in image_paths ]

    def calc_overlap(self):
        """
        座標リストから隣の画像の座標との差を計算
        """
        from statistics import mode
        # x座標リスト
        x_list = sorted([lh.tile_info["x"] for lh in self.labelholders])

        # 隣接する要素間の差を計算して新しいリストを生成
        differences = [x_list[i+1] - x_list[i] for i in range(len(x_list)-1)]
        differences = [ gap for gap in differences if gap > 0]
        # 最頻値取得
        self.gap = mode(differences)
        self.overlap = self.gap // self.downsample
        
    def register_croppedCoord(self):
        """
        Center cropした場合の座標情報を登録
        """
        for lh in self.labelholders:
            lh.tile_info_c = croppedCoord(lh.tile_info, self.overlap)
            lh.tile_info_c = downsampleCoord(lh.tile_info_c)
            
    def prepare_zeroTile(self, num_classes):
        """
        中心画像と周囲画像8枚から成る画像サイズ分のゼロアレイ作成
        さらにそこに貼り付ける場合の各画像の貼り付け座標も記録
        """
        trim = max(2, self.overlap//32)
        
        # 周囲画像の画像開始座標を算出 trimされている版
        start1 = int(0)
        start1_t = trim
        start2 = int(self.tile_size-self.overlap)
        start2_t = start2 + trim
        start3 = int(self.tile_size*2 - self.overlap*2)
        start3_t = start3 + trim
    
        # 各要素の組み合わせ作成
        start_set = [(y_start,x_start) for x_start in [start1_t,start2_t,start3_t] for y_start in [start1_t,start2_t,start3_t]]
        # 番号と開始座標の辞書
        start_dict = dict(zip(range(9), start_set))
        
        # 中心画像の開始座標
        start_dict["center"] = start2
        
        self.trim = trim
        
        # ゼロアレイ作成
        tile = np.zeros(shape=(start3+self.tile_size,start3+self.tile_size,num_classes,9)) 
        return tile, start_dict
    
    def register_surrond(self):
        """
        1つのlabelholderに周囲8枚のlabelholderの情報を紐づける機能
        """
        # defaultdictを使うとkeyが無いときにエラーが出ない
        from collections import defaultdict
    
        # xy座標とimageholderの対応情報を事前に作っておく
        coord_dict = defaultdict(list)
        for lh in self.labelholders:
            x, y = lh.tile_info["x"], lh.tile_info["y"]
            coord_dict[(x, y)].append(lh)
        
        # 1つのlabelholderずつ処理
        for lh in self.labelholders:
            # 中心画像の座標
            x, y = lh.tile_info["x"], lh.tile_info["y"]   

            surround_dict = {}

            # 周囲8枚の画像を座標情報から探す。
            counter = 0
            for x_coord in [x-self.gap, x, x+self.gap]:
                for y_coord in [y-self.gap, y, y+self.gap]:
                    # 中心画像の時はskip
                    if counter != 4:
                        # ImageHolderが入ったリストからその座標のImageHolderを取り出す。返り値はリスト。
                        surround = coord_dict[(x_coord,y_coord)]

                        # 周囲画像がある場合
                        if surround:
                            surround_dict[counter] = surround[0]
                            surround[0].total_call_num += 1
                    counter += 1

            lh.surround = surround_dict
    
    def process(self, model, batch_size, dataset=None,
                centering=True, average_probability=True):
        """
        以下をバッチサイズ単位で逐次的に実行する
        1. 分類モデルで予測
        2. 予測ラベルと周辺予測ラベルの重複領域の確率を平均し、採用ラベルを決定
        3. whole slideサイズのゼロアレイへの貼り付け
        """
        print(f"----- {self.slidename} -----")


        print(average_probability)
        if average_probability:
            if self.num_classes is None:
                self.num_classes = model.get_layer(index=-1).output.shape[-1]
            zeroTile, start_dict = self.prepare_zeroTile(self.num_classes)
            self.register_surrond()
        
        if dataset is None:
            dataset = GeneratePredictDataset_labelholder(
                self.labelholders, batch_size, self.tile_size, centering
            )
        
        complete_list = []
        uncomplete_list = []

        # whole slideへの貼り付けに関するプログレスバー
        pbar = tqdm(total=len(self.labelholders))

        idx = 0    
        for images in tqdm(dataset):
            # バッチサイズで予測処理
            labels = model.predict_on_batch(images)

            # labelholderにlabelを登録
            for label in labels:
                if average_probability:
                    self.labelholders[idx].label = label
                    # overlap/32だけ外縁をトリミング
                    self.labelholders[idx].cropLabel(self.trim)                     
                    # 未処理リストへ追加
                    uncomplete_list.append(self.labelholders[idx])
                else:
                    self.labelholders[idx].label_final = np.argmax(label, axis=-1)
                    # 処理完了リストへ
                    complete_list.append(self.labelholders[idx])
                idx += 1
            
            if average_probability:
                ## 1 batchが終わったらタイル確率平均化->ラベル決定
                for lh in uncomplete_list:
                    if lh.label_final is None:
                        # 周囲の画像のラベルが作成済みかチェック
                        if lh.checkSurround():
                            # 周囲の画像と重複領域の確率を平均 -> ラベル採用
                            lh.averageProb(start_dict, self.num_classes, zeroTile)

                # 必要数呼び出された画像はcomplete_listに移動
                i = 0; end = len(uncomplete_list)
                while i < end:
                    lh = uncomplete_list.pop(0)    

                    if lh.call_num == lh.total_call_num and lh.label_final is not None:
                        lh.removeLabel()
                        complete_list.append(lh)
                    else:
                        uncomplete_list.append(lh)                
                    i+=1

            # complete_listに入ったものをwhole slideに貼り付け
            while complete_list:
                lh = complete_list.pop()
                # Center crop
                lh.centerCrop(self.overlap)
                label = lh.label_final              
                x,y = lh.tile_info_c["x"], lh.tile_info_c["y"]
                self.slide[y:y+label.shape[0],x:x+label.shape[1]] = label
                pbar.update(1)
        pbar.close()
        
        # 問題なければlabelholdersを削除してメモリ解放
        if len(uncomplete_list) > 0:
            print(f"Warning. {len(uncomplete_list)} still remain")
        else:
            del self.labelholders
    
        return self.slide

    def predict_onetile(self, model, 
                        path,
                        centering=True,
                        average_probability=True):
        """
        画像1枚だけをチェックする機能
        """
        # labelholdersから検索
        hit = [ lh for lh in self.labelholders if lh.path == path]
        if not hit:
            tile_info = getTileCoord(path)
            hit = [ lh for lh in self.labelholders if lh.tile_info == tile_info]
        if not hit:
            print("Tile is not found.")
            return

        if average_probability:
            if self.num_classes is None:
                self.num_classes = model.get_layer(index=-1).output.shape[-1]
                
            zeroTile, start_dict = self.prepare_zeroTile(self.num_classes)
            self.register_surrond()            
            surround = list(hit[0].surround.values())
            hit.extend(surround)

        dataset = GeneratePredictDataset_labelholder(
                hit, len(hit), self.tile_size, centering
            )
        labels = model.predict(dataset)
        for lh,label in zip(hit,labels):
            lh.label = label

        if average_probability:
            # overlap/32だけ外縁をトリミング
            [lh.cropLabel(self.trim) for lh in hit]
            hit[0].averageProb(start_dict, self.num_classes, zeroTile)
            for lh in hit:
                lh.call_num = 0
            return hit[0].label_final
        else:
            return np.argmax(hit[0].label,axis=-1)
    




def registerSurrond(labelholders):
    # defaultdictを使うとkeyが無いときにエラーが出ない
    from collections import defaultdict
    
    gap, overlap = calcOverlap(labelholders=labelholders)
    
    # xy座標とimageholderの対応情報を事前に作っておく
    coord_dict = defaultdict(list)
    for lh in labelholders:
        x, y = lh.tile_info["x"], lh.tile_info["y"]
        coord_dict[(x, y)].append(lh)
    
    for lh in tqdm(labelholders):
        # 中心画像の座標
        x, y = lh.tile_info["x"], lh.tile_info["y"]   
    
        surround_dict = {}
        
        # 周囲8枚の画像を座標情報から探す。
        counter = 0
        for x_coord in [x-gap, x, x+gap]:
            for y_coord in [y-gap, y, y+gap]:
                # 中心画像の時はskip
                if counter != 4:
                    # ImageHolderが入ったリストからその座標のImageHolderを取り出す。返り値はリスト。
                    surround = coord_dict[(x_coord,y_coord)]

                    # 周囲画像がある場合
                    if surround:
                        surround_dict[counter] = surround[0]
                        surround[0].total_call_num += 1
                counter += 1
                
        lh.surround = surround_dict

def registerRegion(labelholders, region_size=50000):
    """
    whole slideでの大まかなregionをlabelholderに登録
    デフォルトは5万px * 5万pxの領域
    """
    for labelholder in tqdm(labelholders):
        region_x = labelholder.tile_info["x"] // region_size
        region_y = labelholder.tile_info["y"] // region_size
        labelholder.region = f"{region_x}_{region_y}"
        
# averageProb機能に必要な内容を用意する機能
def prepare_aveProb(
    labelholders,
    trim,
    num_classes,
    region_size=50000,
    image_dir=None,
    extension="jpg",
    image_paths=None,  
):
    """
    Argments
        trim: 予測ラベルの外縁を何ピクセルトリミングするか。Noneならoverlap//32
    
    Return:
        start_dict: trim後のタイル開始座標。centerはtrim前のタイル開始座標
        tile: 9枚の画像を貼り付ける先のゼロアレイ。trim前の形状で作成
    """
    if labelholders is None:
        if image_paths is None:
            image_paths = prepare_image_paths(image_dir, extension)
        labelholders = [ LabelHolder(p) for p in image_paths]
    
    gap, overlap = calcOverlap(labelholders=labelholders)
    tile_size = labelholders[0].tile_info_d["w"]
    
    # 各タイルに周囲タイル情報を登録
    if labelholders[0].surround is None:
        registerSurrond(labelholders)
    # 各タイルのwholeslideでの領域を登録
    registerRegion(labelholders, region_size)
    
    if trim is None:
        trim = overlap // 32
        
    # 周囲画像の画像開始座標を算出 trimされている版
    start1 = int(0); start1_t = trim
    start2 = int(tile_size - overlap); start2_t = start2 + trim
    start3 = int(tile_size*2 - overlap*2); start3_t = start3 + trim
    
    # 各要素の組み合わせ作成
    start_set = [(y_start,x_start) for y_start in [start1_t,start2_t,start3_t] for x_start in [start1_t,start2_t,start3_t]]
    # 番号と開始座標の辞書
    start_dict = dict(zip(range(9), start_set))
    
    # 中心画像のトリム前の開始座標
    start_dict["center"] = start2
    
    # ゼロアレイ作成
    tile = np.zeros(shape=(start3+tile_size,start3+tile_size,num_classes,9))  
    
    return labelholders, start_dict, trim, overlap, tile_size, num_classes, tile        

def AverageProbmap(model,
                   labelholders,
                   batch_size,
                   num_classes,
                   centering=True,
                   trim=None,
                   verbose=True,
                  ):
    """
    labelholderのリストを
    予測 -> 重複領域の確率を平均 -> 最終ラベル化
    する機能
    
    Arguments
        model: 予測に使用するモデル
        labelholders: 同じwhole slide由来のlabelholderのリスト
        batch_size: 予測処理時のバッチ数
        num_classes: 分類クラス数
        centering: 入力明視野画像の値を-1~1にするかどうか。Falseなら0~1に補正。
        trim: 予測ラベルのトリム幅。上下左右からトリミングする。Noneならoverlap//32の値を使用。
        verbose: Trueなら各繰り返し処理のprogress barを表示
    
    Return
        最終ラベルが決定されたlabelholderのリスト。
    """
    
    # 必要なinputの用意、labelholdersのメンバ変数の更新（registerSurround, registerRegion）
    labelholders, start_dict, trim, overlap, tile_size, num_classes, tile = prepare_aveProb(labelholders,trim,num_classes)
    
    # Region情報取得
    regions = sorted(set([ lh.region for lh in labelholders]))
    
    # 結果の保存先
    complete_list = []
    uncomplete_list = []

    # region単位で繰り返し
    for region in tqdm(regions, disable= not verbose):
        # regionのlabelholderのみを抽出
        _labelholders = [ lh for lh in labelholders if lh.region == region]
        # Dataset作成
        dataset = GeneratePredictDataset_labelholder(
            _labelholders, batch_size, tile_size, centering
        )
        
        idx = 0
        for images in tqdm(dataset, disable= not verbose):
            # バッチサイズで予測処理
            labels = model.predict_on_batch(images)
            
            # labelholderにlabelを登録
            for label in labels:
                _labelholders[idx].label = label
                # crop
                _labelholders[idx].cropLabel(trim)
                
                idx += 1

        ## 1 regionが終わったらタイル確率平均化->ラベル決定
        # uncomplete_listに追加
        uncomplete_list.extend(_labelholders)
        for lh in tqdm(uncomplete_list, disable= not verbose):
            if lh.label_final is None:
                # 周囲の画像のラベルが作成済みかチェック
                if lh.checkSurround():
                    # 周囲の画像と重複領域の確率を平均 -> ラベル採用
                    lh.averageProb(start_dict, num_classes, tile)
        
        # 必要数呼び出された画像はcomplete_listに移動
        i = 0; end = len(uncomplete_list)
        while i < end:
            # uncomplete_list先頭のlabelholderを取り出し
            lh = uncomplete_list.pop(0)    
            
            if lh.call_num == lh.total_call_num and lh.label_final is not None:
                # ラベルを削除
                lh.removeLabel()
                # complete_listに追加
                complete_list.append(lh)
            else:
                # uncomplete_list末尾に追加
                uncomplete_list.append(lh)                
            i+=1
        
        print(len(complete_list))
        print(len(uncomplete_list))
    
    return complete_list

def AverageProbmap_wholeslide(
    model,
    labelholders,
    batch_size,
    num_classes,
    centering=True,
    trim=None
):
    # whole slideサイズのゼロアレイ作成
    slide = createEmptySlide(labelholders=labelholders)
    
    # 必要なinputの用意、labelholderのメンバ変数の更新（registerSurround, registerRegion）
    labelholders, start_dict, trim, overlap, tile_size, num_classes, tile = prepare_aveProb(labelholders,trim,num_classes)
    dataset = GeneratePredictDataset_labelholder(
        labelholders, batch_size, tile_size, centering
    )
    
    complete_list = []
    uncomplete_list = []
    
    # whole slideへの貼り付けに関するプログレスバー
    pbar = tqdm(total=len(labelholders))
    
    idx = 0    
    for images in tqdm(dataset):
        # バッチサイズで予測処理
        labels = model.predict_on_batch(images)

        # labelholderにlabelを登録
        for label in labels:
            labelholders[idx].label = label
            # crop
            labelholders[idx].cropLabel(trim)
            # 未処理リストへ追加
            uncomplete_list.append(labelholders[idx])
            idx += 1

        ## 1 batchが終わったらタイル確率平均化->ラベル決定
        for lh in uncomplete_list:
            if lh.label_final is None:
                # 周囲の画像のラベルが作成済みかチェック
                if lh.checkSurround():
                    # 周囲の画像と重複領域の確率を平均 -> ラベル採用
                    lh.averageProb(start_dict, num_classes, tile)

        # 必要数呼び出された画像はcomplete_listに移動
        i = 0; end = len(uncomplete_list)
        while i < end:
            lh = uncomplete_list.pop(0)    

            if lh.call_num == lh.total_call_num and lh.label_final is not None:
                lh.removeLabel()
                complete_list.append(lh)
            else:
                uncomplete_list.append(lh)                
            i+=1
        
        # complete_listに入ったものをwhole slideに貼り付け
        while complete_list:
            lh = complete_list.pop()
            label = lh.label_final
            x,y = lh.tile_info_d["x"], lh.tile_info_d["y"]
            slide[y:y+label.shape[0],x:x+label.shape[0]] = label
            pbar.update(1)
            
    # 問題なければlabelholdersを削除してメモリ解放
    if uncomplete_list:
        print(f"Warning. {len(uncomplete_list)} still remain")
    else:
        del labelholders
    
    return slide