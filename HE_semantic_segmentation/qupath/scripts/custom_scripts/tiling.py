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
    
    differences = [x_list[i+1] - x_list[i] for i in range(len(x_list)-1)]
    differences = [ gap for gap in differences if gap > 0]
    gap = mode(differences)

    downsample = coord_list[0]["d"]
    
    overlap = gap // downsample
    
    return gap, overlap


def getTileCoord(
    path
):
    """
    Extracts tile coordinate information from the file path.
    """
    filename = os.path.basename(path)
    parts = re.split("[,\[\]]", filename)
    
    tile_info = {}
    keys = ["d","x","y","w","h"]
    if not "d=" in path:
        tile_info["d"] = int(1)
        keys.pop(0)
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
    Computes the maximum x and y coordinates from a list of file paths.
    """

    if image_paths is None:
        image_paths = prepare_image_paths(image_dir, extension)
        
    x=[]; y=[]
    for path in image_paths:
        tile_info = getTileCoord(path)
        x.append(tile_info["x"])
        y.append(tile_info["y"])
        
    x_max = max(x); y_max = max(y)
    
    return x_max, y_max


def createEmptySlide(
    image_paths=None,
    labelholders=None,
    downsampled=True
):
    """
    Creates an empty array with whole-slide dimensions.
    """
    if image_paths is None:
        image_paths = [ lh.path for lh in labelholders ]
        
    x_max, y_max = getOriginalImageSize(image_paths=image_paths)
    
    tile_info = getTileCoord(image_paths[0])
    
    tile_size = tile_info["w"]
    d = tile_info["d"] 
    width = math.ceil((x_max+tile_size)/d) if downsampled else x_max+tile_size
    height = math.ceil((y_max+tile_size)/d) if downsampled else y_max+tile_size
    emptySlide = np.zeros(shape=(height,width), dtype=np.uint8)
    
    return emptySlide


def cropLabel(
    label,
    overlap=None,
    trim=None
):
    """
    Crops the tile from all four sides.
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
    Updates coordinate information after cropping all four sides of the tile.
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
    Adjusts coordinate information according to the downsampling factor.
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
    from collections import defaultdict
    
    gap, overlap = calcOverlap(path_list=label_paths)
    
    coord_dict = defaultdict(list)
    for path in label_paths:
        tile_info = getTileCoord(path)
        x, y = tile_info["x"], tile_info["y"]
        coord_dict[(x, y)].append(path)
    
    surround_dict = {}
    
    for center_path in tqdm(label_paths):
        tile_info = getTileCoord(center_path)
        x, y = tile_info["x"], tile_info["y"]   
    
        _surround_dict = {}
        
        counter = 0
        for x_coord in [x-gap, x, x+gap]:
            for y_coord in [y-gap, y, y+gap]:
                surround = coord_dict[(x_coord,y_coord)]
                
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
    start1 = int(0)
    start2 = int(tile_size-overlap)
    start3 = int(tile_size*2-overlap*2)
    
    start_set = [(y_start,x_start) for y_start in [start1,start2,start3] for x_start in [start1,start2,start3]]
    start_dict = dict(zip(range(9), start_set)) 
    tile = np.zeros(shape=(start3+tile_size, start3+tile_size, num_classes,9))
    
    return tile, start_dict


def aveProbTile(
    path,
    tile_size,
    overlap,
    tile,
    surround_dict,
    start_dict,
    num_classes=None,
):   
    if num_classes is None:
        num_classes = np.load(path).shape[-1]      

    trim = overlap//32
    
    _surround_dict = surround_dict[path]
    for counter, surround_path in _surround_dict.items():
        label = np.load(surround_path)
        label = cropLabel(label,trim=trim)
        y_start,x_start = start_dict[counter]
        tile[y_start+trim:y_start+tile_size-trim, x_start+trim:x_start+tile_size-trim,:,counter] = label
        
    center, _ = start_dict[4]
    tile = tile[center:center+tile_size, center:center+tile_size, ...]
    
    tile = np.einsum('ijkl->ijk', tile)/tile.shape[-1]
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
    if label_paths is None:
        label_paths = prepare_image_paths(label_dir,extension)

    if averageProb:
        if num_classes is None:
            num_classes = np.load(label_paths[0]).shape[-1]
            
        surround_dict = prepare_surrond_dict(label_paths)
        emptyTile, start_dict = prepare_emptyTile(
            tile_size, overlap, num_classes
        )

    slide = createEmptySlide(label_paths,downsampled=downsampled)
    
    tile_info_dict = {p: getTileCoord(p) for p in label_paths}

    tile_info_dict = {p: croppedCoord(tile_info,overlap) for p,tile_info in tile_info_dict.items()}
    
    if downsampled:
        tile_info_dict = {p: downsampleCoord(tile_info) for p,tile_info in tile_info_dict.items()}

    if verbose:
        from tqdm.auto import tqdm
        iterator = tqdm(tile_info_dict.items())
    else:
        iterator = tile_info_dict.items()
    for label_path, tile_info in tqdm(tile_info_dict.items()):
        if averageProb:
            tile = aveProbTile(path=label_path,
                               tile_size=tile_size,
                               overlap=overlap,
                               tile=emptyTile,
                               surround_dict=surround_dict,
                               start_dict=start_dict,
                               num_classes=num_classes) 
        else:
            tile = cv2.imread(label_path, cv2.IMREAD_GRAYSCALE)
        
        tile = cropLabel(tile, overlap)
            
        if not downsampled:
            tile = scipy.ndimage.zoom(input=tile, zoom=tile_info["d"], order=0)
        
        x, y = tile_info["x"], tile_info["y"]
        slide[y:y+tile.shape[0],x:x+tile.shape[1]] = tile

    return slide


