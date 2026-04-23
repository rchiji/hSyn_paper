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
    non_background_pixels = np.where(slide>0)
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
    
    min_y, max_y, min_x, max_x = getBoundingBox(slide)
    
    heightSplit = np.arange(min_y, max_y, tile_size)
    heightSplit = np.append(heightSplit , (heightSplit[-1] + (max_y - min_y) % tile_size))
    widthSplit = np.arange(min_x, max_x, tile_size)
    widthSplit = np.append(widthSplit , (widthSplit[-1] + (max_x - min_x) % tile_size))

    for i, w in enumerate(widthSplit[:-1]):
        for j, h in enumerate(heightSplit[:-1]):
            tile = slide[h:heightSplit[j+1], w:widthSplit[i+1]]
            tileHeight = tile.shape[0]
            tileWidth = tile.shape[1]

            io.imsave(f"{save_dir}/{slide_name} [d={downsample},x={w*downsample},y={h*downsample},w={tileWidth*downsample},h={tileHeight*downsample}].png", tile, check_contrast=False)


def write_QuPath_geojson(
    slide,
    classifications,
    downsample,
    slide_name,
    save_dir = "geojson",
    min_area=100,
    x_correction=0,
    y_correction=0,
):
    os.makedirs(save_dir, exist_ok=True)
    
    num_classes = np.max(slide) + 1
    dict_list = []
    for i in tqdm(range(1, num_classes)):
        tmp = np.zeros_like(slide)
        tmp[slide==i] = 1
    
        contours, _ = cv2.findContours(image=tmp, 
                                       mode=cv2.RETR_LIST, 
                                       method=cv2.CHAIN_APPROX_SIMPLE)
    
        contours = [ c for c in contours if c.shape[0] > 2]

        if min_area != 0:
            contours = [contour for contour in contours if cv2.contourArea(contour) >= min_area]
        
        contours = [ np.squeeze(contour).tolist() for contour in contours]
        
        for contour in contours:
            contour.append(contour[0])
    
        contours = [[[c[0]*downsample+x_correction, c[1]*downsample+y_correction] for c in coord] for coord in contours]

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
        
    final_dict = {
        "type":"FeatureCollection",
        "features": dict_list
    }
        
    save_path = f"{save_dir}/{slide_name}.zip"
    with ZipFile(save_path, "w") as z:
        with z.open(f"{slide_name}.geojson", "w") as c:
            c.write(geojson.dumps(final_dict).encode("utf-8"))