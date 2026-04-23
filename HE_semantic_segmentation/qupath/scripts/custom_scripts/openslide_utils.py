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
    Returns tile coordinates for inference from WSI data.

    Arguments:
    - slide
        OpenSlide object of the WSI or file path to the WSI
    - tile_size
        Pixel size of each tile edge used for machine learning
    - overlap
        Overlap width between adjacent tiles
    - downsample_factor
        Downsampling factor applied for machine learning

    Returns:
        List of dictionaries containing coordinate information
    """

    if isinstance(slide,str):
        slide = openslide.open_slide(slide)

    slide_width, slide_height = slide.dimensions

    tile_size_original = tile_size * downsample_factor
    overlap_original = overlap * downsample_factor
    
    start_x=0; start_x_list=[0]
    while start_x + tile_size_original < slide_width:
        start_x = int(start_x  + tile_size_original - overlap_original)
        start_x_list.append(start_x)
        
    start_y=0; start_y_list=[0]
    while start_y + tile_size_original < slide_height:
        start_y = int(start_y  + tile_size_original - overlap_original)
        start_y_list.append(start_y)

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
        with open(json_filepath, "rb") as zip_file:
            zip_data = io.BytesIO(zip_file.read())
        
        with zipfile.ZipFile(zip_data, "r") as zip_ref:
            json_file_name = zip_ref.namelist()[0]
            
            with zip_ref.open(json_file_name) as json_file:
                roi_data = json.load(json_file)
                
    elif "gzip" == mine_encoding:
        with gzip.open(json_filepath, 'rt') as f:
            roi_data = json.load(f)
    else:
        with open(json_filepath, 'r') as json_file:
            roi_data = json.load(json_file)

    return roi_data


def get_ROIList(roi_data):
    if isinstance(roi_data, dict):
        roi_data = roi_data["features"]
        
    ROI_list = []
    for i in range(len(roi_data)):
        coord = roi_data[i]["geometry"]["coordinates"][0]
        ROI_list.append(coord)

    return ROI_list


def check_overlap(
    ROI,
    Tile,
    Threshold = 0,
):
    """
    Determines whether a tile sufficiently overlaps with the region(s) of interest (ROI).
    Arguments:
     - ROI: Shapely object representing the ROI. If multiple, provide as a list.
     - Tile: Shapely object representing the tile
     - Threshold: Minimum area ratio of ROI within the tile (between 0 and 1)
    Returns:
     Boolean indicating whether the tile overlaps the target region with an area ratio above the threshold
    """
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
    Threshold=0,
):
    """
    Extracts only the tiles that contain ROI regions using coordinate information from QuPath Annotation Objects (GeoJSON).
    Arguments:
    - json_filepath
        Path to the GeoJSON file exported from QuPath
    - tile_coords
        List of tile coordinates. Output of the get_tile_coords function.
    - tile_size
        Pixel size of each tile edge used for machine learning
    - downsample_factor
        Downsampling factor applied for machine learning
    - Threshold
        Minimum area ratio of ROI within the tile (between 0 and 1)
    """

    tile_size_original = tile_size * downsample_factor

    roi_data = read_geojson(json_filepath=json_filepath)
    ROI_list = get_ROIList(roi_data)
    ROI_list = [Polygon(coord) for coord in ROI_list]

    final_tile_coords = []
    for tile_coord in tqdm(tile_coords):
        x, y = tile_coord["x"], tile_coord["y"]
        x_end = x + tile_coord["w"]; y_end = y + tile_coord["h"]
        tile_polygon = Polygon([(x,y),(x_end,y),(x_end,y_end),(x,y_end)])
        
        if check_overlap(ROI_list,tile_polygon,Threshold):
            final_tile_coords.append(tile_coord)
        
    return final_tile_coords
