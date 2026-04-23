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
        self.call_num = 0
        self.total_call_num = 0
        self.surround = None
        self.region = None
        self.trimmed = False


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
        self.prepare_labelholders()
        
        super().calc_overlap()
        self.slide = self.create_empty_slide()
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
        x_max = max( tile_info["x"] for tile_info in self.tile_coords)
        y_max = max( tile_info["y"] for tile_info in self.tile_coords) 
        
        width = math.ceil((x_max+self.tile_size_original)/self.downsample)
        height = math.ceil((y_max+self.tile_size_original)/self.downsample)
        
        empty_slide = np.zeros(shape=(height,width), dtype=np.uint8)
        
        return empty_slide

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
        image_batch = [self.read_region(lh.tile_info["x"],lh.tile_info["y"]) for lh in labelholders_batch]
        image_batch = [rescale_image(image, centering) for image in image_batch]
        image_batch = tf.stack(image_batch)
        
        labels = model.predict_on_batch(image_batch)       
        
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

        pbar = tqdm(total=len(self.labelholders))

        idx = 0         
        for i in tqdm(range((len(self.labelholders)-1)//batch_size + 1)):
            start_i = i*batch_size; end_i = start_i + batch_size
            labelholders_batch = self.labelholders[start_i:end_i]

            self.predict_batch(model,labelholders_batch,centering)
            
            for lh in labelholders_batch:
                if average_probability:
                    lh.cropLabel(self.trim)
                    uncomplete_list.append(lh)
                else:
                    lh.label_final = np.argmax(lh.label, axis=-1)
                    complete_list.append(lh)
                
            if average_probability:
                for lh in uncomplete_list:
                    if lh.label_final is None:
                        if lh.checkSurround():
                            lh.averageProb(start_dict, num_classes, zeroTile)

                i = 0; end = len(uncomplete_list)
                while i < end:
                    lh = uncomplete_list.pop(0)    

                    if lh.call_num == lh.total_call_num and lh.label_final is not None:
                        lh.removeLabel()
                        complete_list.append(lh)
                    else:
                        uncomplete_list.append(lh)                
                    i+=1

            while complete_list:
                lh = complete_list.pop()
                lh.centerCrop(self.overlap)
                label = lh.label_final              
                x,y = lh.tile_info_c["x"], lh.tile_info_c["y"]
                self.slide[y:y+label.shape[0],x:x+label.shape[1]] = label
                pbar.update(1)
        pbar.close()
        
        if len(uncomplete_list) > 0:
            print(f"Warning. {len(uncomplete_list)} still remain")
        else:
            del self.labelholders
    
        return self.slide