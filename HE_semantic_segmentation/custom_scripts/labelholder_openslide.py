import os
import numpy as np
from tqdm.auto import tqdm

from datasets import *
from tiling import *
from labelholder import *
from openslide_utils import *


class LabelHolder_OpenSlide(LabelHolder):
    def __init__(self, tile_info):
        self.tile_info = tile_info
        self.tile_info_d = downsampleCoord(self.tile_info)
        self.label = None
        self.label_final = None
        self.call_num = 0 # 何回呼ばれたか記録
        self.total_call_num = 0 # 合計何回呼ばれるかを事前計算した値
        self.surround = None
        self.region = None
        self.trimmed = False # crop済みかの記録  


class LabelHolders_OpenSlide(LabelHolders):
    def __init__(self,
                 slide_path,
                 tile_size=512,
                 overlap=256,
                 downsample_factor=2,
                 json_filepath=None,
                 threshold=0,
                ):
        
        self.slidename = os.path.basename(slide_path)
        self.wsi = openslide.open_slide(slide_path)
        self.tile_size = tile_size
        self.tile_size_original = tile_size * downsample_factor
        self.overlap = overlap
        self.downsample = downsample_factor
        self.json_filepath = json_filepath
        self.threshold = threshold
        # labelholderのリストを登録
        self.prepare_labelholders()
        
        # 隣接タイルとのoverlapや貼り付け座標のずれを計算
        super().calc_overlap()
        # whole slideサイズのゼロアレイ作成
        self.slide = self.create_empty_slide()
        # Center cropした場合の貼り付け先座標の更新
        super().register_croppedCoord()
        
    def prepare_labelholders(self):
        tile_coords = prepare_tile_coords(
            self.wsi, self.tile_size, self.overlap, self.downsample
        )
        if self.json_filepath is not None:
            tile_coords = extract_Tiles(
                self.json_filepath, tile_coords, self.tile_size, self.downsample, self.threshold
            )
        self.tile_coords = tile_coords
        self.labelholders = [LabelHolder_OpenSlide(tile_info) for tile_info in tile_coords]

    def create_empty_slide(self):
        # 最大x, yを調べる
        x_max = max( tile_info["x"] for tile_info in self.tile_coords)
        y_max = max( tile_info["y"] for tile_info in self.tile_coords) 
        
        # whole slideの縦横ピクセル数計算
        width = math.ceil((x_max+self.tile_size_original)/self.downsample)
        height = math.ceil((y_max+self.tile_size_original)/self.downsample)
        
        # ゼロアレイ作成
        empty_slide = np.zeros(shape=(height,width), dtype=np.uint8)
        
        return empty_slide

    # 画像読み込み機能
    def read_region(self,x,y):
        image = self.wsi.read_region(
            location=(x,y),
            level=0,
            size=(self.tile_size_original,self.tile_size_original)
        )
        image = tf.keras.utils.img_to_array(image.convert("RGB"))
        image = tf.convert_to_tensor(image)
        image = tf.image.resize(images=image,
                                size=(self.tile_size,self.tile_size),
                                method=tf.image.ResizeMethod.BILINEAR)
        image = tf.cast(image, tf.float32)        
        return image

    def predict_batch(self,model,labelholders_batch,centering):
        # ImageHolderから明視野画像を取り出してリスト化
        image_batch = [self.read_region(lh.tile_info["x"],lh.tile_info["y"]) for lh in labelholders_batch]
        # 値の範囲変更
        image_batch = [rescale_image(image, centering) for image in image_batch]
        # batch軸で画像を重ねる
        image_batch = tf.stack(image_batch)
        
        # 画像リストを一括で予想
        labels = model.predict_on_batch(image_batch)       
        
        # labelholderにlabelを登録
        with tf.device("/cpu:0"):
            for i in range(len(labels)):
                labelholders_batch[i].label = labels[i]

        del image_batch

    def process(self, model, batch_size, num_classes=None,
                centering=True, average_probability=True):
        
        if average_probability:
            if num_classes is None:
                num_classes = model.get_layer(index=-1).output.shape[-1]
            zeroTile, start_dict = super().prepare_zeroTile(num_classes)
            super().register_surrond()
            
        complete_list = []
        uncomplete_list = []

        # whole slideへの貼り付けに関するプログレスバー
        pbar = tqdm(total=len(self.labelholders))

        idx = 0         
        # batch処理。batch_sizeで割り切れる数に1足すと繰り返し回数が求められる。
        for i in tqdm(range((len(self.labelholders)-1)//batch_size + 1)):
            # 1バッチでの要素範囲
            start_i = i*batch_size; end_i = start_i + batch_size
            # batch size分のImageHolderをリストから取り出す。
            labelholders_batch = self.labelholders[start_i:end_i]

            # 予測とラベル画像の登録
            self.predict_batch(model,labelholders_batch,centering)
            
            for lh in labelholders_batch:
                if average_probability:
                    # overlap/32だけ外縁をトリミング
                    lh.cropLabel(self.trim)
                    # 未処理リストへ追加
                    uncomplete_list.append(lh)
                else:
                    lh.label_final = np.argmax(lh.label, axis=-1)
                    # 処理完了リストへ
                    complete_list.append(lh)
                
            if average_probability:
                ## 1 batchが終わったらタイル確率平均化->ラベル決定
                for lh in uncomplete_list:
                    if lh.label_final is None:
                        # 周囲の画像のラベルが作成済みかチェック
                        if lh.checkSurround():
                            # 周囲の画像と重複領域の確率を平均 -> ラベル採用
                            lh.averageProb(start_dict, num_classes, zeroTile)

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