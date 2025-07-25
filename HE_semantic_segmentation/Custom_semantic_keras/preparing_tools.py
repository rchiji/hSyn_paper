## 221109 update
## 230415 pad_imageを追加。tileの切れ端みたいな明視野画像があった場合に、事前にローカルの画像をpaddingしておく機能。

import os
from tqdm.auto import tqdm
import re
import matplotlib.pyplot as plt
from glob import glob
import numpy as np
import cv2


def MaskToIndex(
    mask_dir:str="Labels",
    save_dir:str="Labels_idx",
    num_classes:int=None,
    path_list:list[str]=None,
    extension:str="png",
    save_ColorMap_jpg:bool=True,
    save_ColorMap_npy:bool=True):
    """
    QuPathで書き出したRGBカラーのマスク画像(png)をクラス数のindex値に変換する。
    背景の白色がindex値0になるように指定。
    色とindex値の対応もindex.jpgとindex.npyとして書き出される。
    最大クラス数を指定した方が処理が早く終わる。
    """

    # 保存先フォルダ作成
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    # mask画像へのパス一覧を事前に取得していない場合はmask画像が入っているフォルダパスをmask_dir引数に指定。
    if path_list is None:
        paths = sorted(glob(f"{mask_dir}/*{extension}"))
        print(f"Get mask path from --{mask_dir}--")
    else:
        paths = sorted(path_list)
    
    print("=============================================")
    
    # すべてのマスク画像を一旦読み込む
    print("Loading masks")
    masks = []
    for path in tqdm(paths):
        mask = cv2.imread(path)
        mask = cv2.cvtColor(mask, cv2.COLOR_BGR2RGB)
        masks.append(mask)
    masks = np.array(masks)
    print("Loading masks --> Done!!")
    print("=============================================")
    
    ## MakeColorMap関数でRGBとlabel_indexの対応dictを取得
    col_idx = MakeColorMap(masks=masks,
                           num_classes=num_classes,
                           save_ColorMap_jpg=save_ColorMap_jpg, 
                           save_ColorMap_npy=save_ColorMap_npy)
    print("=============================================")
    
    ## 対応dictを使ってRGBをindex値に変換したものを保存していく。
    print("Converting masks to index")
    for i in tqdm(range(len(masks))):
        tmp = np.dstack((masks[i], np.zeros(masks[i].shape[:2], 'uint8'))).view('uint32').squeeze(-1)
        # dictのkeyとvalueを使ってRGB積をcolor indexに変換
        for key, value in col_idx.items():
            tmp[tmp==key]=value

        # 保存先フォルダに保存
        name = os.path.basename(paths[i])
        cv2.imwrite(f"./{save_dir}/{name}", tmp.astype("uint8"))
    print(f"Indexed masks were saved to --{save_dir}-- !!")
    print("=============================================")


def MakeColorMap(
    masks:list,
    mask_dir:str=None,
    path_list:list[str]=None,
    extension:str="png",
    num_classes:int=None,
    save_ColorMap_jpg:bool=True,
    save_ColorMap_npy:bool=True):
    """
    ColorMap (マスクのRGBとlabel_indexの対応表)を作成
    返り値はdict型オブジェクト。ColorMapを可視化用に表示する他、jpgやnpyで保存可能。
    白pxがlabel_indexの0 (輝度0 px)に割り当てられるようにしている。
    uniqueなRGB値の探索は時間がかかるので、num_classes引数で最大クラス数を指定しておけば、その数だけの色が見つかり次第処理が終了する。
    """

    # mask画像一覧を事前に取得していない場合はフォルダパスもしくはマスク画像へのパス一覧から画像を読み込む。
    if masks is None and path_list is not None:
        paths = sorted(path_list)
    elif masks is None and mask_dir is not None:
        paths = sorted(glob(f"{mask_dir}/*{extension}"))
        print(f"Getting mask path from {mask_dir}")
        print("=============================================")

    if masks is None and "paths" in locals(): # <-- ローカル変数にpathsがあれば画像を読み込む
        # すべてのマスク画像を一旦読み込む
        print("Loading masks")
        masks = []
        for path in tqdm(paths):
            mask = cv2.imread(path)
            mask = cv2.cvtColor(mask, cv2.COLOR_BGR2RGB)
            masks.append(mask)
        masks = np.array(masks)
        print("Loading masks --> Done!!")
        print("=============================================")
    
    ## ユニークなRGBパターンを見つける。
    print("Getting unique color")
    ch_uniq = []
    if num_classes is None:
        for i in tqdm(range(len(masks))):
            # 2次元画像を1次元化して、uniqueなRGB値のみをch_uniqリストに追記していく。
            ch_uniq.extend(np.unique(masks[i].reshape(-1, masks[i].shape[2]), axis=0))
        # 全画像を通してのuniqueなRGBのみにする。
        ch_uniq = np.unique(ch_uniq, axis=0)
    
    # num_classes引数の指定があれば、それに達した時点でユニークなRGB値の探索は終了する。
    else:
        for i in tqdm(range(len(masks))):
            # 2次元画像を1次元化して、uniqueなRGB値のみをch_uniqリストに追記していく。
            ch_uniq.extend(np.unique(masks[i].reshape(-1, masks[i].shape[2]), axis=0))
            # 全画像を通してのuniqueなRGBのみにする。
            curr_numClass = np.unique(ch_uniq, axis=0)
            if len(curr_numClass) == num_classes:
                ch_uniq = curr_numClass
                break
        
    print(f"{len(ch_uniq)} class were found")
    print("=============================================")

    col_idx = np.expand_dims(ch_uniq, axis = 0)
    # RGB値からuniqueなスカラ値みたいなのを出す。RGBの積っぽい値だけど0があっても0にはならない値。
    col_idx = np.dstack((col_idx, np.zeros(col_idx.shape[:2], 'uint8'))).view("uint32").squeeze(-1)
    col_idx = col_idx.squeeze().tolist() # list型に戻す。
    order = [i[0] for i in sorted(enumerate(col_idx), reverse=True, key=lambda x:x[1])] # 白pxが先頭に来るように並び替えた時のindex順を保存しておく。
    col_idx = sorted(col_idx, reverse=True) # 白pxが先頭に来るように並び替え

    col_idx = dict(zip(col_idx, range(len(col_idx))))
    
    ## ColorMapを可視化して、jpgとして保存
    ch_uniq = ch_uniq[order] # 並び替え
    img = np.expand_dims(ch_uniq, axis = 0)
    plt.figure(figsize=(10,10))
    for i in range(img.shape[1]):
        plt.subplot(1, img.shape[1], i+1)
        plt.title(i)
        plt.axis("off")
        plt.imshow(np.expand_dims(img[:,i,:], axis=0))
    if save_ColorMap_jpg==True:
        plt.savefig("index.jpg")
        print("Color map is saved to current directory!!")
        print("=============================================")
    
    if save_ColorMap_npy==True:
        np.save("index", ch_uniq)
    
    return col_idx

def pad_image(image_paths=None, image_dir=None, extension = "jpg", tileSize = 512):
    if image_paths is None:
        image_paths = glob(f"{image_dir}/*{extension}")

    count = 0
    for path in tqdm(image_paths):
        image = plt.imread(path)
        if image.shape[0] < tileSize or image.shape[1] < tileSize:
            new_image = np.zeros(shape = (tileSize, tileSize, 3), dtype = np.uint8)
            new_image[0:image.shape[0],0:image.shape[1], :] = image
            plt.imsave(path, new_image)
            count += 1
    print(f"{count} images were updated")