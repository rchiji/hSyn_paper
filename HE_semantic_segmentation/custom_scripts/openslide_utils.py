OPENSLIDE_PATH = r"C:\Users\admin\Downloads\software\openslide-win64-20231011\bin"

import os
if hasattr(os, 'add_dll_directory'):
    # Windows
    with os.add_dll_directory(OPENSLIDE_PATH):
        import openslide
else:
    import openslide
    
import openslide
import mimetypes
import io
import zipfile
import json
import gzip
from shapely.geometry import Polygon
from tqdm.auto import tqdm


# 画像読み込み機能
def read_region(slide,location,size):
    image = slide.read_region(location=location,
                              level=0,
                              size=size)
    image = tf.keras.utils.img_to_array(image.convert("RGB"))
    image = tf.convert_to_tensor(image)
    image = tf.cast(image, tf.float32)
    
    return image 
    

def prepare_tile_coords(
    slide,
    tile_size:int = 512,
    overlap:int = 256,
    downsample_factor:int = 2
):
    """
    WSIデータから予測用タイル画像の座標を取得する機能

    Arguments:
    - slide
        WSIのOpenSlideオブジェクト または WSIのファイルパス
    - tile_size
        機械学習に使用するタイル辺のピクセルサイズ
    - overlap
        タイルの重複幅
    - downsample_factor
        機械学習に使用する際のダウンサンプリング係数
    
    Return:
     座標情報の辞書のリスト
    """

    # slide引数がファイルパスならOpenSlideオブジェクトを作成
    if isinstance(slide,str):
        slide = openslide.open_slide(slide)

    # WSLデータのピクセルサイズを取得
    slide_width, slide_height = slide.dimensions

    # 元の解像度の場合のピクセルサイズ
    tile_size_original = tile_size * downsample_factor
    overlap_original = overlap * downsample_factor
    
    # タイルの開始座標を算出。縦方向、横方向それぞれで行う。
    # 元のスライドサイズを超えるまでタイルの開始座標を記録していく。
    start_x=0; start_x_list=[0]
    while start_x + tile_size_original < slide_width:
        start_x = int(start_x  + tile_size_original - overlap_original)
        start_x_list.append(start_x)
        
    start_y=0; start_y_list=[0]
    while start_y + tile_size_original < slide_height:
        start_y = int(start_y  + tile_size_original - overlap_original)
        start_y_list.append(start_y)

    # タイル開始座標を[x,y]の組み合わせたリストに。さらに全座標を一つのリストに。
    tile_coords = []
    for i in range(len(start_x_list)):
        for j in range(len(start_y_list)):
            tile_info = {"d":downsample_factor,
                         "x":start_x_list[i],
                         "y":start_y_list[j],
                         "w":tile_size_original,
                         "h":tile_size_original}
            tile_coords.append(tile_info)

    return tile_coords


def read_geojson(json_filepath):
    mime_type, mine_encoding = mimetypes.guess_type(url = json_filepath)
    if mime_type is not None and "zip" in mime_type:
        # ZIPファイルをメモリに読み込む
        with open(json_filepath, "rb") as zip_file:
            zip_data = io.BytesIO(zip_file.read())
        
        # ZIPファイルを解凍せずにJSONデータを直接読み込む
        with zipfile.ZipFile(zip_data, "r") as zip_ref:
            # ZIPファイル内のJSONファイル名
            json_file_name = zip_ref.namelist()[0]
            
            with zip_ref.open(json_file_name) as json_file:
                # JSONファイルを読み込む
                roi_data = json.load(json_file)
                
    elif "gzip" == mine_encoding:
        with gzip.open(json_filepath, 'rt') as f:
            # JSONデータを読み込む
            roi_data = json.load(f)
    else:
        with open(json_filepath, 'r') as json_file:
            roi_data = json.load(json_file)

    return roi_data


def get_ROIList(roi_data):
    if isinstance(roi_data, dict):
        roi_data = roi_data["features"]
        
    ROI_list = []
    ## ROIの数だけ繰り返し処理。
    for i in range(len(roi_data)):
        # 座標情報を取り出す
        coord = roi_data[i]["geometry"]["coordinates"][0]
        ROI_list.append(coord)

    return ROI_list


def check_overlap(
    ROI,
    Tile,
    Threshold = 0,
):
    """
    Tile内にROIが含まれるか判断する機能
    
    Arguments:
     - ROI: ROIのsharpyオブジェクト。複数の場合はリストで渡す。
     - Tile: Tileのsharpyオブジェクト
     - Threshold: Tile中のROIの割合の閾値。0-1の間で指定。

    Return:
     予測タイル座標が予測対象エリアに閾値以上の面積割合で跨るかどうかのbool値
    """
    # ROI引数の指定がROIのリストの場合
    if isinstance(ROI, list):
        checklist = []
        for roi in ROI:
            area = Tile.intersection(roi).area
            ratio = area/Tile.area
            checklist.append(ratio > Threshold)
            if any(checklist) is True:
                break                
        return any(checklist)            
        
    else:
        area = Tile.intersection(roi).area
        ratio = area/Tile.area        
        return ratio > Threshold 


def extract_Tiles(
    json_filepath:str,
    tile_coords:list,
    tile_size=512,
    downsample_factor=2,
    Threshold=0, # タイル中のROI割合閾値
):
    """
    QuPathのAnnotation Objectの座標情報（geojson）を使って、ROIを含むタイルだけを抽出する機能。

    Arguments:
    - json_filepath
        QuPathから書き出したgeojsonファイルのパス
    - tile_coords
        タイル画像の座標リスト。get_tile_coords機能の出力。
    - tile_size
        機械学習に使用するタイル辺のピクセルサイズ
    - downsample_factor
        機械学習に使用する際のダウンサンプリング係数
    - Threshold
        Tile中のROIの割合の閾値。0-1の間で指定。
    """

    # ダウンサンプリング前のWSI座標でのタイルサイズ
    tile_size_original = tile_size * downsample_factor

    # geojsonファイルからsharpy Polygonオブジェクトを作成
    roi_data = read_geojson(json_filepath=json_filepath)
    ROI_list = get_ROIList(roi_data)
    ROI_list = [Polygon(coord) for coord in ROI_list]

    # 最終的なタイル座標の保存先リスト
    final_tile_coords = []
    # 各タイル座標のPolygonオブジェクトがQuPath ROIを含むかどうか確認
    for tile_coord in tqdm(tile_coords):
        x, y = tile_coord["x"], tile_coord["y"]
        x_end = x + tile_coord["w"]; y_end = y + tile_coord["h"]
        tile_polygon = Polygon([(x,y),(x_end,y),(x_end,y_end),(x,y_end)])
        
        if check_overlap(ROI_list,tile_polygon,Threshold):
            final_tile_coords.append(tile_coord)
        
    return final_tile_coords
