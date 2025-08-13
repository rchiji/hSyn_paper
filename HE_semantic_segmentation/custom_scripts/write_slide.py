import os
import numpy as np
from glob import glob
from tqdm.auto import tqdm
import scipy.ndimage
from skimage import io
import json
import geojson
from zipfile import ZipFile

from datasets import *
from tiling import *
from labelholder import *


def getBoundingBox(
    slide
):
    # 背景以外の領域を抽出
    non_background_pixels = np.where(slide>0)
    # BoundingBoxの座標
    min_y, min_x = np.min(non_background_pixels, axis=1)
    max_y, max_x = np.max(non_background_pixels, axis=1)

    return min_y, max_y, min_x, max_x


def writeLabelSlide(
    slide,
    slide_name,
    downsample,
    save_dir,
):
    os.makedirs(save_dir, exist_ok=True)
    
    min_y, max_y, min_x, max_x = getBoundingBox(slide)
    slide = slide[min_y:max_y, min_x:max_x]
    
    x = min_x * downsample
    y = min_y * downsample
    width =  max_x * downsample - x
    height = max_y * downsample - y
    save_name = f"{slide_name} [d={downsample},x={x},y={y},w={width},y={height}]"

    io.imsave(f"{save_dir}/{save_name}.png", slide, check_contrast=False)


def downsampleSlide(
    slide,
    downsample_factor:int
):
    slide = scipy.ndimage.zoom(slide, 1/downsample_factor, order=0)
    return slide


def saveSplitedLabel(
    slide,
    tile_size,
    downsample,
    slide_name,
    save_dir,
):
    os.makedirs(save_dir,exist_ok=True)
    
    # 画像サイズを取得
    min_y, max_y, min_x, max_x = getBoundingBox(slide)
    
    # 分割時の座標リスト作成
    heightSplit = np.arange(min_y, max_y, tile_size)
    heightSplit = np.append(heightSplit , (heightSplit[-1] + (max_y - min_y) % tile_size)) # 剰余はappendで追加
    widthSplit = np.arange(min_x, max_x, tile_size)
    widthSplit = np.append(widthSplit , (widthSplit[-1] + (max_x - min_x) % tile_size))

    for i, w in enumerate(widthSplit[:-1]):
        for j, h in enumerate(heightSplit[:-1]):
            tile = slide[h:heightSplit[j+1], w:widthSplit[i+1]]
            tileHeight = tile.shape[0]
            tileWidth = tile.shape[1]

            # QuPathの書き出し時の名前っぽく保存
            io.imsave(f"{save_dir}/{slide_name} [d={downsample},x={w*downsample},y={h*downsample},w={tileWidth*downsample},h={tileHeight*downsample}].png", tile, check_contrast=False)


def write_QuPath_geojson(
    slide,
    classifications, # 背景以外の分類ラベル名のリスト
    downsample,
    slide_name,
    save_dir = "geojson",
    min_area=100, # Objectとして残す最小ピクセル面積
    x_correction=0, # x座標の補正値
    y_correction=0, # y座標の補正値
):
    os.makedirs(save_dir, exist_ok=True)
    
    num_classes = np.max(slide) + 1
    dict_list = []
    for i in tqdm(range(1, num_classes)):
        # 同じ形状のゼロアレイを用意
        tmp = np.zeros_like(slide)
        # 0,1に変換
        tmp[slide==i] = 1
    
        # 輪郭座標の取得
        contours, _ = cv2.findContours(image=tmp, 
                                       mode=cv2.RETR_LIST, 
                                       method=cv2.CHAIN_APPROX_SIMPLE)
    
        # 座標が3点以上のものだけにする
        contours = [ c for c in contours if c.shape[0] > 2]

        # 小さいObjectは除く
        if min_area != 0:
            contours = [contour for contour in contours if cv2.contourArea(contour) >= min_area]
        
        # squeezeして無駄次元を無くす ＆ numpyをlistに変換
        contours = [ np.squeeze(contour).tolist() for contour in contours]
        
        # 輪郭座標の始点を末尾に追加（Objectを閉じるため）
        for contour in contours:
            contour.append(contour[0])
    
        # ダウンサンプル前の座標に更新
        contours = [[[c[0]*downsample+x_correction, c[1]*downsample+y_correction] for c in coord] for coord in contours]

        # QuPathスタイルのObject情報dict作成
        for contour in contours:
            dict_data = {"type":"Feature",
                         "geometry":{"type":"Polygon",
                                     "coordinates":[contour]},
                         "properties":{"objectType":"annotation",
                                       "classification":{
                                           "name":f"{classifications[i-1]}",
                                           "color":[]
                                       },
                                       "isLocked":"false",
                                       "measurements":[]}
                        }
        
            dict_list.append(dict_data)
        
    # 全てのObjectをまとめたdict作成
    final_dict = {
        "type":"FeatureCollection",
        "features": dict_list
    }
        
    # 書き出し
    save_path = f"{save_dir}/{slide_name}.zip"
    with ZipFile(save_path, "w") as z:
        with z.open(f"{slide_name}.geojson", "w") as c:
            c.write(geojson.dumps(final_dict).encode("utf-8"))