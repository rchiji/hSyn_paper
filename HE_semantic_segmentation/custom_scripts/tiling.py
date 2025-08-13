import os
import numpy as np
import re
from glob import glob
from tqdm.auto import tqdm
import scipy.ndimage
import math
import cv2
import copy
import matplotlib.pyplot as plt

try:
    from .datasets import prepare_image_paths
except ImportError:
    from datasets import prepare_image_paths

def calcOverlap(
    path_list=None,
    coord_list=None,
    labelholders=None,    
):
    from statistics import mode
    if coord_list is None and labelholders is None:
        coord_list = [ getTileCoord(p) for p in path_list]
    else:
        coord_list = [ lh.tile_info for lh in labelholders]
        
    x_list = [coord["x"] for coord in coord_list]
    
    # 隣接する要素間の差を計算して新しいリストを生成
    differences = [x_list[i+1] - x_list[i] for i in range(len(x_list)-1)]
    differences = [ gap for gap in differences if gap > 0]
    # 最頻値取得
    gap = mode(differences)

    downsample = coord_list[0]["d"]
    
    overlap = gap // downsample
    
    return gap, overlap


def getTileCoord(
    path
):
    """
    ファイルパスからタイル座標を取得
    """
    filename = os.path.basename(path)
    # , . [ ]でsplitしてリスト化
    parts = re.split("[,\[\]]", filename)
    
    tile_info = {}
    keys = ["d","x","y","w","h"]
    # ファイル名にダウンサンプリング係数が無い場合はd=1にする。
    if not "d=" in path:
        tile_info["d"] = int(1)
        keys.pop(0)
    # ファイル名から各種情報を取得
    for key in keys:
        tmp = [ part for part in parts if part.startswith(f"{key}=")][0]
        tile_info[key] = int(re.split("[=]", tmp)[1])

    return tile_info


def getOriginalImageSize(
    image_paths=None,
    image_dir=None,
    extension=None
):
    """
    ファイル名リストから最大x,y座標を取得
    """

    if image_paths is None:
        image_paths = prepare_image_paths(image_dir, extension)
        
    x=[]; y=[]
    for path in image_paths:
        tile_info = getTileCoord(path)
        x.append(tile_info["x"])
        y.append(tile_info["y"])
        
    # 最大値を取得
    x_max = max(x); y_max = max(y)
    
    return x_max, y_max


def createEmptySlide(
    image_paths=None,
    labelholders=None,
    downsampled=True
):
    """
    whole slideサイズの空アレイ作成
    """
    # 最大x, yを調べる
    if image_paths is None:
        image_paths = [ lh.path for lh in labelholders ]
        
    x_max, y_max = getOriginalImageSize(image_paths=image_paths)
    
    tile_info = getTileCoord(image_paths[0])
    
    # whole slideの縦横ピクセル数計算
    tile_size = tile_info["w"]
    d = tile_info["d"] 
    width = math.ceil((x_max+tile_size)/d) if downsampled else x_max+tile_size
    height = math.ceil((y_max+tile_size)/d) if downsampled else y_max+tile_size
    # ゼロアレイ作成
    emptySlide = np.zeros(shape=(height,width), dtype=np.uint8)
    
    return emptySlide


def cropLabel(
    label,
    overlap=None,
    trim=None
):
    """
    タイルの上下左右をトリミングする機能
    """
    if trim is None:
        trim = overlap // 2
    label = label[trim:label.shape[0]-trim, trim:label.shape[1]-trim]
    return label


def croppedCoord(
    tile_info,
    overlap
):
    """
    座標情報をタイルの上下左右をトリミングした後の座標に更新する機能
    """
    downsample = tile_info["d"]
    overlap = overlap * downsample
    trim = overlap // 2

    new_tile_info = {
        "d": tile_info["d"],
        "x": tile_info["x"] + trim,
        "y": tile_info["y"] + trim,
        "w": tile_info["w"] - trim*2,
        "h": tile_info["h"] - trim*2
    }
    return new_tile_info


def downsampleCoord(
    tile_info,
):
    """
    ダウンサンプリング係数で座標情報を補正する機能
    """
    downsample = tile_info["d"]
    new_tile_info = {
        "d": tile_info["d"],
        "x": tile_info["x"] // downsample,
        "y": tile_info["y"] // downsample,
        "w": tile_info["w"] // downsample,
        "h": tile_info["h"] // downsample
    }

    return new_tile_info


def prepare_surrond_dict(label_paths):
    # defaultdictを使うとkeyが無いときにエラーが出ない
    from collections import defaultdict
    
    gap, overlap = calcOverlap(path_list=label_paths)
    
    # xy座標とimageholderの対応情報を事前に作っておく
    coord_dict = defaultdict(list)
    for path in label_paths:
        tile_info = getTileCoord(path)
        x, y = tile_info["x"], tile_info["y"]
        coord_dict[(x, y)].append(path)
    
    surround_dict = {}
    
    for center_path in tqdm(label_paths):
        # 中心画像の座標
        tile_info = getTileCoord(center_path)
        x, y = tile_info["x"], tile_info["y"]   
    
        _surround_dict = {}
        
        # 周囲8枚の画像を座標情報から探す。中心画像も含める。
        counter = 0
        for x_coord in [x-gap, x, x+gap]:
            for y_coord in [y-gap, y, y+gap]:
                # ImageHolderが入ったリストからその座標のImageHolderを取り出す。返り値はリスト。
                surround = coord_dict[(x_coord,y_coord)]
                
                # 周囲画像がある場合
                if surround:
                    _surround_dict[counter] = surround[0]
            counter += 1
                
        surround_dict[center_path] = _surround_dict
    
    return surround_dict


def prepare_emptyTile(
    tile_size,
    overlap,
    num_classes
):
    # 9枚の画像を貼り合わせ先を作成
    # 周囲画像の画像開始座標を算出
    start1 = int(0)
    start2 = int(tile_size-overlap)
    start3 = int(tile_size*2-overlap*2)
    
    # 各要素の組み合わせ作成
    start_set = [(y_start,x_start) for y_start in [start1,start2,start3] for x_start in [start1,start2,start3]]
    # 番号と開始座標の辞書
    start_dict = dict(zip(range(9), start_set)) 
    # ゼロアレイ作成
    tile = np.zeros(shape=(start3+tile_size, start3+tile_size, num_classes,9))
    
    return tile, start_dict


def aveProbTile(
    path, # 中心画像
    tile_size, # QuPathで指定したtile px数
    overlap, # QuPathで指定したoverlap px数
    tile,
    surround_dict,
    start_dict,
    num_classes=None,
):   
    if num_classes is None:
        num_classes = np.load(path).shape[-1]      

    # トリミング幅
    trim = overlap//32
    
    # 周囲画像ファイルパスと貼り付け位置の辞書取り出し
    _surround_dict = surround_dict[path]
    for counter, surround_path in _surround_dict.items():
        label = np.load(surround_path)
        # overlap pxの1/32を上下左右からトリミング
        label = cropLabel(label,trim=trim)
        y_start,x_start = start_dict[counter]
        tile[y_start+trim:y_start+tile_size-trim, x_start+trim:x_start+tile_size-trim,:,counter] = label
        
    # 中心画像の領域をcrop
    center, _ = start_dict[4]
    tile = tile[center:center+tile_size, center:center+tile_size, ...]
    
    # 周囲画像の画像とで確率を平均 (tile_size,tile_size,num_classes,9) -> (tile,tile,num_classes)
    tile = np.einsum('ijkl->ijk', tile)/tile.shape[-1]
    # 最大確率のラベルを採用 (tile_size,tile_size,num_classes) -> (tile_size,tile_size)
    tile = np.argmax(tile,axis=-1)

    return tile


def tilingLabel(
    label_dir=None,
    overlap=256,
    tile_size=512,
    averageProb:bool=True,
    downsampled=True,
    extension=None,
    label_paths=None,
    num_classes=None,
):
    # ラベル画像のパス一覧を取得
    if label_paths is None:
        label_paths = prepare_image_paths(label_dir,extension)

    if averageProb:
        if num_classes is None:
            num_classes = np.load(label_paths[0]).shape[-1]
            
        # 周囲画像ファイルパスと貼り付け位置の辞書
        surround_dict = prepare_surrond_dict(label_paths)
        # 確率を平均する際の一時貼り付け先
        emptyTile, start_dict = prepare_emptyTile(
            tile_size, overlap, num_classes
        )

    ## 貼り付け先のゼロアレイ作成
    slide = createEmptySlide(label_paths,downsampled=downsampled)
    
    # ラベル画像パスと座標情報の辞書作成
    tile_info_dict = {p: getTileCoord(p) for p in label_paths}

    # Center Cropした場合の座標情報に更新
    tile_info_dict = {p: croppedCoord(tile_info,overlap) for p,tile_info in tile_info_dict.items()}
    
    # ダウンサンプリングされたwhole slideの座標情報に更新
    if downsampled:
        tile_info_dict = {p: downsampleCoord(tile_info) for p,tile_info in tile_info_dict.items()}

    # ラベルの読み込みと貼り付け
    if verbose:
        from tqdm.auto import tqdm
        iterator = tqdm(tile_info_dict.items())
    else:
        iterator = tile_info_dict.items()
    for label_path, tile_info in tqdm(tile_info_dict.items()):
        if averageProb:
            # 重複領域の確率平均 -> ラベル採用
            tile = aveProbTile(path=label_path,
                               tile_size=tile_size,
                               overlap=overlap,
                               tile=emptyTile,
                               surround_dict=surround_dict,
                               start_dict=start_dict,
                               num_classes=num_classes) 
        else:
            # ラベル画像の読み込み
            tile = cv2.imread(label_path, cv2.IMREAD_GRAYSCALE)
        
        # Center crop
        tile = cropLabel(tile, overlap)
            
        # オリジナルサイズで行う場合は、タイルをアップサンプリング
        if not downsampled:
            tile = scipy.ndimage.zoom(input=tile, zoom=tile_info["d"], order=0)
        
        # 貼り付け先の座標
        x, y = tile_info["x"], tile_info["y"]
        # 貼り付け
        slide[y:y+tile.shape[0],x:x+tile.shape[1]] = tile

    return slide


