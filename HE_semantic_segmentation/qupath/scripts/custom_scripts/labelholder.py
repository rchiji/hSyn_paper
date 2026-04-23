import os
import numpy as np
from tqdm.auto import tqdm
import scipy.ndimage
import cv2

from datasets import *
from tiling import *
from prediction import *

class LabelHolder:
    def __init__(self, path):
        self.path = path
        self.tile_info = getTileCoord(self.path)
        self.tile_info_d = downsampleCoord(self.tile_info)
        self.label = None
        self.label_final = None
        self.call_num = 0
        self.total_call_num = 0
        self.surround = None
        self.region = None
        self.trimmed = False
    
    def checkSurround(self):
        """
        Checks whether prediction labels for all surrounding tiles are available.
        """
        check = [surround.label is not None for surround in self.surround.values()]
        if check:
            return all(check)   
        else:
            return False
    
    def cropLabel(self,trim=None, overlap=None):
        """
        Crops the label image by the specified trim width from all four sides.
        """
        if trim is None:
            trim = overlap // 32
        if self.trimmed is False:
            self.label = self.label[trim:self.label.shape[0]-trim, trim:self.label.shape[1]-trim]
            self.trimmed = True
    
    def centerCrop(self, overlap):
        """
        Crops the final label image by half of the overlap width from all four sides.
        """
        trim = overlap // 2
        self.label_final = self.label_final[trim:self.label_final.shape[0]-trim, trim:self.label_final.shape[1]-trim]
        

    def averageProb(self, start_dict, num_classes, tile):
        """
        Averages label probabilities across the overlapping regions of the center tile and its surrounding tiles.
        """
        tile_size = self.tile_info_d["w"]
        
        y_start, x_start = start_dict[4]
        tile[y_start:y_start+self.label.shape[0], x_start:x_start+self.label.shape[1],:,4] = self.label
        
        for counter, slh in self.surround.items():
            y_start, x_start = start_dict[counter] 
            label = slh.label
            tile[y_start:y_start+label.shape[0], x_start:x_start+label.shape[1],:,counter] = label
            slh.call_num +=1
        
        center = start_dict["center"]
        tile = tile[center:center+tile_size, center:center+tile_size]
        
        tile = np.einsum('ijkl->ijk', tile)/num_classes
        tile = np.argmax(tile,axis=-1)
    
        self.label_final = tile
        
    def removeLabel(self):
        """
        Deletes the label image if it has already been referenced the required number of times.
        """
        if self.call_num == self.total_call_num:
            del self.label


class LabelHolders:
    def __init__(self,
                 labelholders=None,                  
                 image_paths=None,
                 image_dir=None,
                ):
        
        self.labelholders = labelholders
        self.image_paths = image_paths
        self.image_dir = image_dir
        if labelholders is None:           
            self.prepare_labelholders()
        
        self.filename = os.path.basename(self.labelholders[0].path)
        self.slidename = re.split("[ ]", self.filename)[0]
        self.downsample = self.labelholders[0].tile_info["d"]
        self.tile_size = self.labelholders[0].tile_info_d["w"]
        self.num_classes = None
        self.calc_overlap()
        self.slide = createEmptySlide(labelholders=self.labelholders)
        self.register_croppedCoord()

    
    def prepare_labelholders(self):
        """
        Creates a list of LabelHolder instances from either a list of brightfield tile image paths or a directory containing the images.
        """
        if self.image_paths is None:
            image_paths = prepare_image_paths(self.image_dir)
            
        image_paths = sorted(image_paths)
        self.labelholders = [ LabelHolder(path) for path in image_paths ]

    def calc_overlap(self):
        """
        Computes the spatial gap between adjacent tiles from coordinate information.
        """
        from statistics import mode
        x_list = sorted([lh.tile_info["x"] for lh in self.labelholders])

        differences = [x_list[i+1] - x_list[i] for i in range(len(x_list)-1)]
        differences = [ gap for gap in differences if gap > 0]
        self.gap = mode(differences)
        self.overlap = self.gap // self.downsample
        
    def register_croppedCoord(self):
        """
        Registers coordinate information after center cropping.
        """
        for lh in self.labelholders:
            lh.tile_info_c = croppedCoord(lh.tile_info, self.overlap)
            lh.tile_info_c = downsampleCoord(lh.tile_info_c)
            
    def prepare_zeroTile(self, num_classes):
        """
        Creates a zero-initialized array for a 9-tile layout (center tile + 8 neighbors), and records the placement coordinates for each tile.
        """
        trim = max(2, self.overlap//32)
        
        start1 = int(0)
        start1_t = trim
        start2 = int(self.tile_size-self.overlap)
        start2_t = start2 + trim
        start3 = int(self.tile_size*2 - self.overlap*2)
        start3_t = start3 + trim
    
        start_set = [(y_start,x_start) for x_start in [start1_t,start2_t,start3_t] for y_start in [start1_t,start2_t,start3_t]]
        start_dict = dict(zip(range(9), start_set))
        
        start_dict["center"] = start2
        
        self.trim = trim
        
        tile = np.zeros(shape=(start3+self.tile_size,start3+self.tile_size,num_classes,9)) 
        return tile, start_dict
    
    def register_surrond(self):
        """
        Associates each LabelHolder with its 8 neighboring LabelHolder instances.
        """
        from collections import defaultdict
    
        coord_dict = defaultdict(list)
        for lh in self.labelholders:
            x, y = lh.tile_info["x"], lh.tile_info["y"]
            coord_dict[(x, y)].append(lh)
        
        for lh in self.labelholders:
            x, y = lh.tile_info["x"], lh.tile_info["y"]   

            surround_dict = {}

            counter = 0
            for x_coord in [x-self.gap, x, x+self.gap]:
                for y_coord in [y-self.gap, y, y+self.gap]:
                    if counter != 4:
                        surround = coord_dict[(x_coord,y_coord)]

                        if surround:
                            surround_dict[counter] = surround[0]
                            surround[0].total_call_num += 1
                    counter += 1

            lh.surround = surround_dict
    
    def process(self, model, batch_size, dataset=None,
                centering=True, average_probability=True):
        """
        Sequentially performs the following steps in batches:
        1. Run model inference
        2. Average probabilities in overlapping regions between the center tile and surrounding tiles, and determine the final label
        3. Paste the result into a whole-slide-sized array
        """
        print(f"----- {self.slidename} -----")


        print(average_probability)
        if average_probability:
            if self.num_classes is None:
                self.num_classes = model.get_layer(index=-1).output.shape[-1]
            zeroTile, start_dict = self.prepare_zeroTile(self.num_classes)
            self.register_surrond()
        
        if dataset is None:
            dataset = GeneratePredictDataset_labelholder(
                self.labelholders, batch_size, self.tile_size, centering
            )
        
        complete_list = []
        uncomplete_list = []

        pbar = tqdm(total=len(self.labelholders))

        idx = 0    
        for images in tqdm(dataset):
            labels = model.predict_on_batch(images)

            for label in labels:
                if average_probability:
                    self.labelholders[idx].label = label
                    self.labelholders[idx].cropLabel(self.trim)                     
                    uncomplete_list.append(self.labelholders[idx])
                else:
                    self.labelholders[idx].label_final = np.argmax(label, axis=-1)
                    complete_list.append(self.labelholders[idx])
                idx += 1
            
            if average_probability:
                for lh in uncomplete_list:
                    if lh.label_final is None:
                        if lh.checkSurround():
                            lh.averageProb(start_dict, self.num_classes, zeroTile)

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

    def predict_onetile(self, model, 
                        path,
                        centering=True,
                        average_probability=True):
        """
        Performs prediction for a single tile for inspection.
        """
        hit = [ lh for lh in self.labelholders if lh.path == path]
        if not hit:
            tile_info = getTileCoord(path)
            hit = [ lh for lh in self.labelholders if lh.tile_info == tile_info]
        if not hit:
            print("Tile is not found.")
            return

        if average_probability:
            if self.num_classes is None:
                self.num_classes = model.get_layer(index=-1).output.shape[-1]
                
            zeroTile, start_dict = self.prepare_zeroTile(self.num_classes)
            self.register_surrond()            
            surround = list(hit[0].surround.values())
            hit.extend(surround)

        dataset = GeneratePredictDataset_labelholder(
                hit, len(hit), self.tile_size, centering
            )
        labels = model.predict(dataset)
        for lh,label in zip(hit,labels):
            lh.label = label

        if average_probability:
            [lh.cropLabel(self.trim) for lh in hit]
            hit[0].averageProb(start_dict, self.num_classes, zeroTile)
            for lh in hit:
                lh.call_num = 0
            return hit[0].label_final
        else:
            return np.argmax(hit[0].label,axis=-1)
    




def registerSurrond(labelholders):
    from collections import defaultdict
    
    gap, overlap = calcOverlap(labelholders=labelholders)
    
    coord_dict = defaultdict(list)
    for lh in labelholders:
        x, y = lh.tile_info["x"], lh.tile_info["y"]
        coord_dict[(x, y)].append(lh)
    
    for lh in tqdm(labelholders):
        x, y = lh.tile_info["x"], lh.tile_info["y"]   
    
        surround_dict = {}
        
        counter = 0
        for x_coord in [x-gap, x, x+gap]:
            for y_coord in [y-gap, y, y+gap]:
                if counter != 4:
                    surround = coord_dict[(x_coord,y_coord)]

                    if surround:
                        surround_dict[counter] = surround[0]
                        surround[0].total_call_num += 1
                counter += 1
                
        lh.surround = surround_dict

def registerRegion(labelholders, region_size=50000):
    """
    Assigns coarse regions to each LabelHolder on the whole slide.
    The default region size is 50,000 × 50,000 pixels.
    """
    for labelholder in tqdm(labelholders):
        region_x = labelholder.tile_info["x"] // region_size
        region_y = labelholder.tile_info["y"] // region_size
        labelholder.region = f"{region_x}_{region_y}"
        
def prepare_aveProb(
    labelholders,
    trim,
    num_classes,
    region_size=50000,
    image_dir=None,
    extension="jpg",
    image_paths=None,  
):
    """
    Arguments:
        trim: Number of pixels to trim from the border of predicted labels.
              If None, overlap // 32 is used.
    Returns:
        start_dict: Tile start coordinates after trimming (center corresponds to pre-trim coordinates)
        tile: Zero-initialized array for placing 9 tiles (constructed with pre-trim dimensions)
    """
    if labelholders is None:
        if image_paths is None:
            image_paths = prepare_image_paths(image_dir, extension)
        labelholders = [ LabelHolder(p) for p in image_paths]
    
    gap, overlap = calcOverlap(labelholders=labelholders)
    tile_size = labelholders[0].tile_info_d["w"]
    
    if labelholders[0].surround is None:
        registerSurrond(labelholders)
    registerRegion(labelholders, region_size)
    
    if trim is None:
        trim = overlap // 32
        
    start1 = int(0); start1_t = trim
    start2 = int(tile_size - overlap); start2_t = start2 + trim
    start3 = int(tile_size*2 - overlap*2); start3_t = start3 + trim
    
    start_set = [(y_start,x_start) for y_start in [start1_t,start2_t,start3_t] for x_start in [start1_t,start2_t,start3_t]]
    start_dict = dict(zip(range(9), start_set))
    
    start_dict["center"] = start2
    
    tile = np.zeros(shape=(start3+tile_size,start3+tile_size,num_classes,9))  
    
    return labelholders, start_dict, trim, overlap, tile_size, num_classes, tile        

def AverageProbmap(model,
                   labelholders,
                   batch_size,
                   num_classes,
                   centering=True,
                   trim=None,
                   verbose=True,
                  ):
    """
    Processes a list of LabelHolder instances by:
    prediction -> averaging probabilities in overlapping regions -> final label assignment
    Arguments:
        model: Model used for inference
        labelholders: List of LabelHolder instances from the same whole slide
        batch_size: Batch size for inference
        num_classes: Number of classes
        centering: Whether to normalize input images to [-1, 1] (if False, normalized to [0, 1])
        trim: Trim width applied to predicted labels. If None, overlap // 32 is used.
        verbose: If True, displays progress bars
    Returns:
        List of LabelHolder instances with finalized labels.
    """
    
    labelholders, start_dict, trim, overlap, tile_size, num_classes, tile = prepare_aveProb(labelholders,trim,num_classes)
    
    regions = sorted(set([ lh.region for lh in labelholders]))
    
    complete_list = []
    uncomplete_list = []

    for region in tqdm(regions, disable= not verbose):
        _labelholders = [ lh for lh in labelholders if lh.region == region]
        dataset = GeneratePredictDataset_labelholder(
            _labelholders, batch_size, tile_size, centering
        )
        
        idx = 0
        for images in tqdm(dataset, disable= not verbose):
            labels = model.predict_on_batch(images)
            
            for label in labels:
                _labelholders[idx].label = label
                _labelholders[idx].cropLabel(trim)
                
                idx += 1

        uncomplete_list.extend(_labelholders)
        for lh in tqdm(uncomplete_list, disable= not verbose):
            if lh.label_final is None:
                if lh.checkSurround():
                    lh.averageProb(start_dict, num_classes, tile)
        
        i = 0; end = len(uncomplete_list)
        while i < end:
            lh = uncomplete_list.pop(0)    
            
            if lh.call_num == lh.total_call_num and lh.label_final is not None:
                lh.removeLabel()
                complete_list.append(lh)
            else:
                uncomplete_list.append(lh)                
            i+=1
        
        print(len(complete_list))
        print(len(uncomplete_list))
    
    return complete_list

def AverageProbmap_wholeslide(
    model,
    labelholders,
    batch_size,
    num_classes,
    centering=True,
    trim=None
):
    slide = createEmptySlide(labelholders=labelholders)
    
    labelholders, start_dict, trim, overlap, tile_size, num_classes, tile = prepare_aveProb(labelholders,trim,num_classes)
    dataset = GeneratePredictDataset_labelholder(
        labelholders, batch_size, tile_size, centering
    )
    
    complete_list = []
    uncomplete_list = []
    
    pbar = tqdm(total=len(labelholders))
    
    idx = 0    
    for images in tqdm(dataset):
        labels = model.predict_on_batch(images)

        for label in labels:
            labelholders[idx].label = label
            labelholders[idx].cropLabel(trim)
            uncomplete_list.append(labelholders[idx])
            idx += 1

        for lh in uncomplete_list:
            if lh.label_final is None:
                if lh.checkSurround():
                    lh.averageProb(start_dict, num_classes, tile)

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
            label = lh.label_final
            x,y = lh.tile_info_d["x"], lh.tile_info_d["y"]
            slide[y:y+label.shape[0],x:x+label.shape[0]] = label
            pbar.update(1)
            
    if uncomplete_list:
        print(f"Warning. {len(uncomplete_list)} still remain")
    else:
        del labelholders
    
    return slide