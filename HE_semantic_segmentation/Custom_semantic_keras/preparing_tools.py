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
    Converts RGB mask images (PNG) exported from QuPath into class index values.
    The background white color is assigned to index 0.
    The correspondence between colors and index values is also saved as index.jpg and index.npy.
    Specifying the maximum number of classes can speed up processing.
    """

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    if path_list is None:
        paths = sorted(glob(f"{mask_dir}/*{extension}"))
        print(f"Get mask path from --{mask_dir}--")
    else:
        paths = sorted(path_list)
    
    print("=============================================")
    
    print("Loading masks")
    masks = []
    for path in tqdm(paths):
        mask = cv2.imread(path)
        mask = cv2.cvtColor(mask, cv2.COLOR_BGR2RGB)
        masks.append(mask)
    masks = np.array(masks)
    print("Loading masks --> Done!!")
    print("=============================================")
    
    col_idx = MakeColorMap(masks=masks,
                           num_classes=num_classes,
                           save_ColorMap_jpg=save_ColorMap_jpg, 
                           save_ColorMap_npy=save_ColorMap_npy)
    print("=============================================")
    
    print("Converting masks to index")
    for i in tqdm(range(len(masks))):
        tmp = np.dstack((masks[i], np.zeros(masks[i].shape[:2], 'uint8'))).view('uint32').squeeze(-1)
        for key, value in col_idx.items():
            tmp[tmp==key]=value

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
    Creates a ColorMap (mapping table between mask RGB values and label indices).
    Returns a dictionary object. The ColorMap can be displayed for visualization and also saved as JPG or NPY.
    White pixels are assigned to label index 0 (intensity value 0).
    Since searching for unique RGB values can be time-consuming, specifying the maximum number of classes with the num_classes argument allows the process to stop as soon as that number of colors is found.
    """

    if masks is None and path_list is not None:
        paths = sorted(path_list)
    elif masks is None and mask_dir is not None:
        paths = sorted(glob(f"{mask_dir}/*{extension}"))
        print(f"Getting mask path from {mask_dir}")
        print("=============================================")

    if masks is None and "paths" in locals():
        print("Loading masks")
        masks = []
        for path in tqdm(paths):
            mask = cv2.imread(path)
            mask = cv2.cvtColor(mask, cv2.COLOR_BGR2RGB)
            masks.append(mask)
        masks = np.array(masks)
        print("Loading masks --> Done!!")
        print("=============================================")
    
    print("Getting unique color")
    ch_uniq = []
    if num_classes is None:
        for i in tqdm(range(len(masks))):
            ch_uniq.extend(np.unique(masks[i].reshape(-1, masks[i].shape[2]), axis=0))
        ch_uniq = np.unique(ch_uniq, axis=0)
    
    else:
        for i in tqdm(range(len(masks))):
            ch_uniq.extend(np.unique(masks[i].reshape(-1, masks[i].shape[2]), axis=0))
            curr_numClass = np.unique(ch_uniq, axis=0)
            if len(curr_numClass) == num_classes:
                ch_uniq = curr_numClass
                break
        
    print(f"{len(ch_uniq)} class were found")
    print("=============================================")

    col_idx = np.expand_dims(ch_uniq, axis = 0)
    col_idx = np.dstack((col_idx, np.zeros(col_idx.shape[:2], 'uint8'))).view("uint32").squeeze(-1)
    col_idx = col_idx.squeeze().tolist()
    order = [i[0] for i in sorted(enumerate(col_idx), reverse=True, key=lambda x:x[1])]
    col_idx = sorted(col_idx, reverse=True)

    col_idx = dict(zip(col_idx, range(len(col_idx))))
    
    ch_uniq = ch_uniq[order]
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