## 221223 updated
"""
predict_on_batchを一旦全て画像で行うと、予測結果をCPUメモリに保持させることになり、PCスペック次第でメモリ上限まで達してしまう。
1000枚の画像で64GBメモリは足りないぐらい。
そこからCPUメモリの圧縮が起こるので止まりはしないが急に遅くなる。

予測バッチが終わるごとに確率マップ平均化の準備ができた画像を探して、確率マップを平均化する。9回呼ばれると予測マスクをImageHolder属性から消せるため、CPUメモリに溜まらないようにできる。
"""

import os
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'
import cv2
import re
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from tqdm.auto import tqdm
import tensorflow as tf
from tensorflow import keras
from keras import layers
from keras import backend as K
from Custom_semantic_keras.tiling_tools import *
from Custom_semantic_keras.losses import *

custom_objects={'UpdatedMeanIoU': UpdatedMeanIoU,
                "Dice_SCCE_loss": Dice_SCCE_loss,
                "Dice":Dice,
                "IoU":IoU}

def GetTrainedModel(model_path, custom_object=True):
    """""
    訓練済みモデル読み込み機能
    ・model_pathに訓練済みモデル.h5ファイルを指定
    ・UpdatedMeanIoUクラスをモデル訓練時に使用している場合はIoUオプションをTrueに。
    ・h5ファイルからモデルをloadするとカスタムオブジェクトを指定する必要あり。
    https://www.tensorflow.org/tutorials/keras/save_and_load?hl=ja
    """""
    
    if custom_object is True:
        model = keras.models.load_model(model_path, custom_objects=custom_objects)
    else:
        model = keras.models.load_model(model_path)
    return model


def read_pred_image(image_path:str, 
                    imageSize:int):
    """
    画像読み込み機能定義
    """
    image = tf.io.read_file(image_path)
    image = tf.image.decode_png(image)
    image = image[:,:,:3]
    image.set_shape([None, None, 3])
    image = tf.image.resize(images=image, size=[imageSize, imageSize])
    image = image / 127.5 - 1
    return image


def GetSlideNames(path_list:list[str],
                  file_dir:str=None,
                  prefix:str=None,
                  extension:str="jpg"):
    """
    スライド名リスト作成機能
    ・画像パス一覧か画像が入ったディレクトリへのパス
    ・ディレクトリから読み込む場合は指定した拡張子の画像のみを扱う。
    ・prefixとして指定した文字は削除してスライド名として取り出せる。
    """
    # 画像のパス一覧の指定がNoneの場合、dir=引数の指定先から画像のパス一覧を取得する。
    if path_list is None:
        path_list = sorted(glob(f"{file_dir}/*{extension}"))
        print(f"Got mask path from {file_dir}")
    else:
        path_list = sorted(path_list)
    
    print("--- Getting slide names ---")
    slidenames = []
    for path in path_list:
        path = os.path.basename(path) # ファイル名を取得
        if prefix is not None:
            if path.startswith(prefix):
                path = path.replace(prefix, "")
                slidenames.append(re.split("[ ]", path)[0])            
        else:
            slidenames.append(re.split("[ ]", path)[0]) # 半角スペースで区切ってスライド名と座標情報を分ける。
        slidenames = sorted(list(set(slidenames)))
    print(f"{len(slidenames)} slides were found")
    return slidenames

    
## 画像の座標情報などを保存するために、ファイルパスを引数にしたインスタンスを作成
class ImageHolder:
    """""
    画像情報を保持させるクラス
    画像へのパスを引数に指定。
    メンバ変数
    ・画像パス、ファイル名
    ・座標情報、サイズ、ダウンサンプリング係数
    ・パスから読み込んだ画像
    ・訓練済みモデルで予測したマスク画像
    ・確率マップ平均化した後のマスク画像
    """""  
    def __init__(self, path):
        self.path = path
        self.filename = os.path.basename(path)
        self.slidename = re.split("[ ]", self.filename)[0] # 半角スペースで区切ってスライド名を保存
        self.tile_info = self.getTileInfoDict()
        self.tile_size = self.tile_info["w"] # ダウンサンプル前の元の画像pxサイズを記録
        self.image = None
        self.read_pred_image() # <-- init内に画像読み込みを入れておく。
        self.imageSize = self.image.shape[0]
        self.mask = None
        self.mask_final = None
        self.call_num = 0 # 何回呼ばれたか記録
    
    # タイル情報取得機能
    def getTileInfoDict(self):
        # ファイル名から, [], (), 半角スペースで文字列を分割したリスト作成
        parts = [p for p in re.split("[ ,\\[\\]\\(\\)]", self.filename)]
        tile_info = dict()
        # 特定の座標情報を抜きだす。
        tile_info["x"] = int([p for p in parts if p.startswith("x")][0][2:])
        tile_info["y"] = int([p for p in parts if p.startswith("y")][0][2:])
        tile_info["w"] = int([p for p in parts if p.startswith("w")][0][2:])
        tile_info["h"] = int([p for p in parts if p.startswith("h")][0][2:])
        # ダウンサンプリング係数が1の場合はd=がない場合があるので、その対応
        if not "d=" in self.filename:
            tile_info["d"] = int(1)
        else:
            tile_info["d"] = int([p for p in parts if p.startswith("d")][0][2:])
        return tile_info
        
    # 画像読み込み機能
    def read_pred_image(self, center_zero=True):
        """""
        cpuメモリに保持させるためにwith tf.device("/cpu:0")内で画像を読み込む
        """"" 
        if self.image is None:
            with tf.device("/cpu:0"):
                image = tf.io.read_file(self.path)
                image = tf.image.decode_png(image)
                image = image[:,:,:3]
                image.set_shape([None, None, 3])
                image = tf.cast(image, tf.float32)
                image = image/127.5 - 1 if center_zero is True else image/255
                self.image = image
        return self.image 
    
    # その画像が既に9回呼び出されていたら画像データとマスクデータのメンバ変数を消す。
    def check_call_count(self):
        if self.call_num == 9:
            self.image = None
            self.mask = None
        return self # ImageHolderインスタンス自体を返す
    
## 重複領域の確率マップを平均化する機能
def AverageProbmap(model,
                   path_list:list[str],
                   numClass:int=14, 
                   imageSize:int=512, 
                   overlap:int=256, 
                   downsample:int=5,
                   batch_size:int=10,
                   save_tilemask:bool=True,
                   save_dir:str="Pred_masks",
                   model_path:str=None,
                   custom_object:bool=True,
                   softmax:bool=True
                   ):
    """""
    訓練済みモデルを読み込んで、明視野画像から予測マスクを作る一連の機能
    1) 画像パス一覧からImageHolderインスタンス一覧作成
    2) バッチ単位で訓練済みモデルで予測
    3) 重複領域の確立マップを平均化
    4) ローカルにマスク画像を保存 (オプション)
    """""    
    
    ## モデルを読み込んでいなかったら、model_pathからモデル読み込み
    if model is None:
        model = GetTrainedModel(model_path=model_path, custom_object=custom_object)

    ## 保存先フォルダの作成
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
        
    # タイル画像の座標位置ずれ計算 (ダウンサンプリングされる前の値)
    gap = imageSize*downsample - overlap*downsample
        
    ## 失敗した画像パスの一覧をglobal変数として残るようにする。
    global fail_list; fail_list=[]
    
    def checkSurround(imgholder, gap):
        """
        周囲画像に予測後マスクがあるか判定
        """
        # 中心画像の座標
        x = imgholder.tile_info["x"]; y = imgholder.tile_info["y"]
        check_list = []
        
        # 中心を含めた9枚の画像を座標情報から探す。
        for x_coord in [x-gap, x, x+gap]:
            for y_coord in [y-gap, y, y+gap]:
                # ImageHolderが入ったリストからその座標のImageHolderを取り出す。リスト内包表記なので返り値がリストのまま。
                tmp = [im for im in imgholder_list if im.tile_info["x"]==x_coord and im.tile_info["y"]==y_coord]
                # tmpの中身が空ではない場合 (周囲画像がある場合)
                if len(tmp)>0:
                    # 周囲画像の予測マスクがまだ無ければloopを終了
                    if tmp[0].mask is None:
                        # check_listにはFalseと残す
                        check_list.append(tmp[0].mask is not None)
                        break
                    else:
                        # 周囲画像の予測マスクがあればcheck_listにTrueと残す。
                        check_list.append(tmp[0].mask is not None)

            # nested for loopの時にouter loopごと終了させる書き方。
            # 周囲画像はあるが、どれか1枚でも予測マスクが未だ保持してなければ処理を終える。
            else:
                continue
            break
        
        # 残ったcheck_listを判定。
        # all()は空リストでもTrueと出るので、まず空リストかどうか判定してから進める。
        if check_list:
            # リスト内全てがTrueかどうかを返す
            return all(check_list)         
        
    def returnMask(imgholder, imageSize=imageSize, numClass=numClass):
        """
        ImageHolderインスタンスの有無を判定してあればマスク予想して結果を返す。
        なければ背景アレイを返す。
        """
        if imgholder:
            imgholder[0].call_num += 1 # 呼び出し回数を記録
            return imgholder[0].mask
        else:
            dummy = tf.Variable(tf.zeros(shape=(imageSize,imageSize,numClass), dtype=K.floatx()))
            dummy = dummy[...,0].assign(1) # 周囲画像が無い = 組織の淵。背景クラスのはずなので、背景クラス(index 0)に背景の確率値1を入れておく。
            return dummy
        
    def calcAveProb(imgholder, imageSize, overlap, numClass, gap, fail_list):
        """
        画像サイズ*クラス数*9枚の4軸アレイを作成。
        9枚の画像をタイルに配置していく。最後の軸で各画像の次元を分ける。
        """
        # 周囲画像の画像開始座標を算出
        Start1 = int(0)
        Start2 = int(imageSize-overlap)
        Start3 = int(imageSize*2-overlap*2)
    
        try:
            tile = tf.Variable(tf.zeros(shape=(Start3+imageSize, Start3+imageSize, numClass, 9), dtype=K.floatx()))
            
            # 周囲8枚の画像を座標から検索。returnMask関数でImageHolderインスタンスからマスク画像を返す。assignでタイルに配置。
            # imgholder引数を中心画像として扱う。
            x = imgholder.tile_info["x"]; y = imgholder.tile_info["y"]
            counter = 0
            # top, center, bottomの順
            for y_coord, y_start in zip([y-gap, y, y+gap], [Start1, Start2, Start3]):
                # left, center, rightの順
                for x_coord, x_start in zip([x-gap, x, x+gap], [Start1, Start2, Start3]):
                    tile = tile[y_start:y_start+imageSize, x_start:x_start+imageSize,:,counter].assign(returnMask([im for im in imgholder_list if im.tile_info["x"]==x_coord and im.tile_info["y"]==y_coord]))
                    counter += 1

            # 最後の軸で確率マップの平均を取り、argmaxで最大確率のクラスを採用
            tile = tf.argmax(tf.reduce_mean(tile, axis=-1), axis=-1)

            with tf.device("/cpu:0"):
                tile = tile.numpy().astype("int8")
                # 元の中心画像領域をcropして、mask_final属性に保存
                imgholder.mask_final = tile[Start2:Start2+imageSize,Start2:Start2+imageSize]
        
        # エラーがあったらfail_listに画像パスを追記して次の画像へ進む
        except Exception as e:
            print(f"{imgholder.filename} is fail")
            print("== error info =="); print(f"type: {str(type(e))}")
            # print(f"args: {str(e.args)}")
            print(f"error: {str(e)}")
            fail_list.append(imgholder.path)
    
    ##########
    ## 画像パスからImageHolderインスタンス作成。
    imgholder_list = []
    print("--- Registering image info ---")
    for i in tqdm(range(len(path_list))):
        try:
            imgholder_list.append(ImageHolder(path=path_list[i]))
            
        except Exception as e:
            print(f"{path_list[i]} is fail")
            print("== error info =="); print(f"type: {str(type(e))}"); print(f"error: {str(e)}")
            fail_list.append(path_list[i])
            
    ##########
    print(f"--- Predicting images (batch size: {batch_size})---")
    # batch処理。batch_sizeで割り切れる数に1足すと繰り返し回数が求められる。
    for i in tqdm(range((len(imgholder_list)-1)//batch_size + 1)):
        # 1バッチでの要素範囲
        start_i = i*batch_size; end_i = start_i + batch_size
        # batch size分のImageHolderをリストから取り出す。
        imgholder_batch = imgholder_list[start_i:end_i]
        # ImageHolderから明視野画像を取り出してリスト化
        img_batch = [imgholder.image for imgholder in imgholder_batch]
        # batch軸で画像を重ねる
        img_batch = tf.stack(img_batch)
        # 画像リストを一括で予想
        masks = model.predict_on_batch(img_batch) 
        if softmax is True:
            masks = K.softmax(masks, axis=-1).numpy()
            
        with tf.device("/cpu:0"):
            # 予測マスクを元のImageHolderのmask属性に保存
            for i in range(len(masks)):
                imgholder_batch[i].mask = masks[i]
                imgholder_batch[i].image = None # メモリ解放のため明視野画像は削除
        
        #####
        # バッチ単位の予測が終わる度に確率マップの平均化を進める。
        for imgholder in imgholder_list:
            if imgholder.mask_final is None:
                # 周囲の画像の予測マスクが作られてるかチェック
                if checkSurround(imgholder, gap) is True:
                    # 準備できていれば確率マップの平均化を実行
                    calcAveProb(imgholder, imageSize, overlap, numClass, gap, fail_list)
                    
        # 画像が9回呼ばれたimageholderの予測マスク情報を削除
        for im in imgholder_list:
            im.check_call_count()
        
    ##########
    if save_tilemask is True:
        print("--- Saving masks ---")
        with tf.device("/cpu:0"):
            for i in tqdm(range(len(imgholder_list))):
                try:
                    saveTile = imgholder_list[i].mask_final
                    saveName = "Labels_" + imgholder_list[i].filename # ファイル名の先頭にLabels_をつける
                    saveName = saveName[:-4] # 拡張子のところ4文字分を消す
                    saveName = saveName + ".png"
                    cv2.imwrite(f"{save_dir}/{saveName}", saveTile)    
                except:
                    pass

    if len(fail_list)>0:
        print(f"{len(fail_list)} images were fail")
        fail_list = "\n".join(fail_list)
        with open('fail_list.txt', 'a') as f:
            f.write(fail_list)
    
    return imgholder_list


## 重複領域の確率マップを平均化する機能
def AverageProbmapSlides(
    model,
    path_list:list[str]=None,
    image_dir:str="tiles_512jpg",
    extension="jpg",
    numClass:int=14,
    imageSize:int=512,
    overlap:int=256,
    downsample:int=5,
    batch_size:int=10,
    save_tilemask:bool=False,
    save_dir:str="Pred_masks",
    model_path:str=None,
    custom_object:bool=True,
    softmax:bool=True,
    makeColorMask:bool=True,
    colormap_path="index.npy"
    ):
    """""
    AverageProbmapをスライド単位で実行し、そのままタイリングまで実行する。
    """""  
    
    ## モデルを読み込んでいなかったら、model_pathからモデル読み込み
    if model is None:
        model = GetTrainedModel(model_path=model_path, custom_object=custom_object)  
        
    if path_list is None:
        # 画像パス一覧を取得
        imgPath = sorted(glob(f"{image_dir}/*{extension}"))
    else:
        imgPath = sorted(path_list)
    
    # 画像が入っているフォルダや画像パスからスライド名一覧を取得
    slidenames = sorted(GetSlideNames(path_list=imgPath,file_dir=None,prefix=None,extension=extension))
    
    global fail_list; fail_list=[] # failした画像パスを保存する用の空リスト
    
    # slideごとにAverageProbmapを実行。どこかでエラーが出たら次のスライドの処理に進む。
    for slide in tqdm(slidenames):
        try:
            print(f"<<<<<< Predicting {slide} >>>>>>")
            paths = [s for s in imgPath if slide in s]
            print(f"{len(paths)} images exist")
            imgholder_list = AverageProbmap(
                model=model,path_list=paths,numClass=numClass, imageSize=imageSize,overlap=overlap,downsample=downsample,
                batch_size=batch_size, save_tilemask=save_tilemask, save_dir=save_dir, softmax=softmax)

            MergeTiles_from_ImageHolder(
                imgholder_list=imgholder_list, overlap_px=overlap, imageSize=imageSize, downsample=downsample,
                save_dir=save_dir, makeColorMask=makeColorMask, save_ColorMask_dir=f"{save_dir}_color", colormap_path=colormap_path)

            print(f"{slide} is done!!\n")
            K.clear_session()
        
        except Exception as e:
            print(f"{slide} is fail")
            print("== error info =="); print(f"type: {str(type(e))}")
            print(f"args: {str(e.args)}"); print(f"error: {str(e)}")
            
def SimplePrediction(path_list, 
                     save_dir="Pred_masks",
                     model_path = "Model_11class_220916.h5",
                     imageSize=256,
                     custom_object=True,
                     save_format_npy=False):
    """""
    訓練済みモデルを読み込んで、明視野画像から予測マスクを作る一連の機能
    重複領域を考慮せずにシンプルにマスク画像を作成する版
    """""
    ## モデル読み込み
    model = GetTrainedModel(model_path=model_path, custom_object=custom_object)
    
    ## 保存先フォルダの作成
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    """ 明視野画像一覧からマスク画像作成"""
    for i in tqdm(range(len(image_paths))):
        path = image_paths[i]
        tmp = read_pred_image(image_path=path, imageSize=imageSize)
        tmp = model.predict_on_batch(np.expand_dims(tmp, axis=0))
        tmp = np.squeeze(tmp)
        tmp = np.argmax(tmp, axis=-1).astype("int8")
        
        # 保存
        saveName = os.path.basename(path)
        saveName = "Labels_" + saveName # ファイル名の先頭にLabels_をつける
        saveName = saveName[:-4] # 拡張子のところ消す
        
        if save_format_npy==True:
            np.save(f"{save_dir}/{saveName}", tmp)
        else:
            saveName = saveName + ".png"
            cv2.imwrite(f"{save_dir}/{saveName}", tmp)

