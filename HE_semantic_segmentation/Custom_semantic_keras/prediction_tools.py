import os
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'
import cv2
import re
import numpy as np
from glob import glob
import matplotlib.pyplot as plt
from tqdm.auto import tqdm
import tensorflow as tf
from tensorflow import keras
from keras import layers
from keras import backend as K
from Custom_semantic_keras.tiling_tools import *
from Custom_semantic_keras.losses import *

custom_objects={'UpdatedMeanIoU': UpdatedMeanIoU,
                "Dice_SCCE_loss": Dice_SCCE_loss,
                "Dice":Dice,
                "IoU":IoU}

def GetTrainedModel(model_path, custom_object=True):
    """
    Loads a trained model.
    - Specify the trained model .h5 file in model_path
    - If the UpdatedMeanIoU class was used during training, set the IoU option to True
    - When loading a model from an .h5 file, custom objects must be specified
    https://www.tensorflow.org/tutorials/keras/save_and_load
    """
    
    if custom_object is True:
        model = keras.models.load_model(model_path, custom_objects=custom_objects)
    else:
        model = keras.models.load_model(model_path)
    return model


def read_pred_image(image_path:str, 
                    imageSize:int):
    """
    Defines an image loading function.
    """
    image = tf.io.read_file(image_path)
    image = tf.image.decode_png(image)
    image = image[:,:,:3]
    image.set_shape([None, None, 3])
    image = tf.image.resize(images=image, size=[imageSize, imageSize])
    image = image / 127.5 - 1
    return image


def GetSlideNames(path_list:list[str],
                  file_dir:str=None,
                  prefix:str=None,
                  extension:str="jpg"):
    """
    Creates a list of slide names.
    - Accepts either a list of image paths or a directory containing images
    - When loading from a directory, only images with the specified extension are used
    - Characters specified as a prefix can be removed before extracting the slide name
    """
    if path_list is None:
        path_list = sorted(glob(f"{file_dir}/*{extension}"))
        print(f"Got mask path from {file_dir}")
    else:
        path_list = sorted(path_list)
    
    print("--- Getting slide names ---")
    slidenames = []
    for path in path_list:
        path = os.path.basename(path)
        if prefix is not None:
            if path.startswith(prefix):
                path = path.replace(prefix, "")
                slidenames.append(re.split("[ ]", path)[0])            
        else:
            slidenames.append(re.split("[ ]", path)[0])
        slidenames = sorted(list(set(slidenames)))
    print(f"{len(slidenames)} slides were found")
    return slidenames

    
class ImageHolder:
    """
    Class for storing image-related information.
    Takes the image path as an argument.

    Member variables:
    - Image path and file name
    - Coordinate information, size, and downsampling factor
    - Image loaded from the path
    - Mask image predicted by the trained model
    - Mask image after probability map averaging
    """
    def __init__(self, path):
        self.path = path
        self.filename = os.path.basename(path)
        self.slidename = re.split("[ ]", self.filename)[0]
        self.tile_info = self.getTileInfoDict()
        self.tile_size = self.tile_info["w"]
        self.image = None
        self.read_pred_image()
        self.imageSize = self.image.shape[0]
        self.mask = None
        self.mask_final = None
        self.call_num = 0
    
    def getTileInfoDict(self):
        parts = [p for p in re.split("[ ,\\[\\]\\(\\)]", self.filename)]
        tile_info = dict()
        tile_info["x"] = int([p for p in parts if p.startswith("x")][0][2:])
        tile_info["y"] = int([p for p in parts if p.startswith("y")][0][2:])
        tile_info["w"] = int([p for p in parts if p.startswith("w")][0][2:])
        tile_info["h"] = int([p for p in parts if p.startswith("h")][0][2:])
        if not "d=" in self.filename:
            tile_info["d"] = int(1)
        else:
            tile_info["d"] = int([p for p in parts if p.startswith("d")][0][2:])
        return tile_info
        
    def read_pred_image(self, center_zero=True):
        """
        Loads images within with tf.device("/cpu:0") to ensure they are stored in CPU memory.
        """
        if self.image is None:
            with tf.device("/cpu:0"):
                image = tf.io.read_file(self.path)
                image = tf.image.decode_png(image)
                image = image[:,:,:3]
                image.set_shape([None, None, 3])
                image = tf.cast(image, tf.float32)
                image = image/127.5 - 1 if center_zero is True else image/255
                self.image = image
        return self.image 
    
    def check_call_count(self):
        if self.call_num == 9:
            self.image = None
            self.mask = None
        return self
    
def AverageProbmap(model,
                   path_list:list[str],
                   numClass:int=14, 
                   imageSize:int=512, 
                   overlap:int=256, 
                   downsample:int=5,
                   batch_size:int=10,
                   save_tilemask:bool=True,
                   save_dir:str="Pred_masks",
                   model_path:str=None,
                   custom_object:bool=True,
                   softmax:bool=True
                   ):
    """
    A series of functions that loads a trained model and generates predicted masks from brightfield images.
    1) Create a list of ImageHolder instances from a list of image paths
    2) Run prediction with the trained model in batches
    3) Average probability maps in overlapping regions
    4) Save mask images locally (optional)
    """ 
    
    if model is None:
        model = GetTrainedModel(model_path=model_path, custom_object=custom_object)

    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
        
    gap = imageSize*downsample - overlap*downsample
        
    global fail_list; fail_list=[]
    
    def checkSurround(imgholder, gap):
        """
        Checks whether predicted masks are available for the surrounding images.
        """
        x = imgholder.tile_info["x"]; y = imgholder.tile_info["y"]
        check_list = []
        
        for x_coord in [x-gap, x, x+gap]:
            for y_coord in [y-gap, y, y+gap]:
                tmp = [im for im in imgholder_list if im.tile_info["x"]==x_coord and im.tile_info["y"]==y_coord]
                if len(tmp)>0:
                    if tmp[0].mask is None:
                        check_list.append(tmp[0].mask is not None)
                        break
                    else:
                        check_list.append(tmp[0].mask is not None)

            else:
                continue
            break
        
        if check_list:
            return all(check_list)         
        
    def returnMask(imgholder, imageSize=imageSize, numClass=numClass):
        """
        Checks whether an ImageHolder instance is available.
        If present, predicts the mask and returns the result.
        Otherwise, returns a background array.
        """
        if imgholder:
            imgholder[0].call_num += 1
            return imgholder[0].mask
        else:
            dummy = tf.Variable(tf.zeros(shape=(imageSize,imageSize,numClass), dtype=K.floatx()))
            dummy = dummy[...,0].assign(1)
            return dummy
        
    def calcAveProb(imgholder, imageSize, overlap, numClass, gap, fail_list):
        """
        Creates a 4D array with dimensions: image size × number of classes × 9 tiles.
        The 9 images are arranged as tiles, with the last axis representing each image separately.
        """
        Start1 = int(0)
        Start2 = int(imageSize-overlap)
        Start3 = int(imageSize*2-overlap*2)
    
        try:
            tile = tf.Variable(tf.zeros(shape=(Start3+imageSize, Start3+imageSize, numClass, 9), dtype=K.floatx()))
            
            x = imgholder.tile_info["x"]; y = imgholder.tile_info["y"]
            counter = 0
            for y_coord, y_start in zip([y-gap, y, y+gap], [Start1, Start2, Start3]):
                for x_coord, x_start in zip([x-gap, x, x+gap], [Start1, Start2, Start3]):
                    tile = tile[y_start:y_start+imageSize, x_start:x_start+imageSize,:,counter].assign(returnMask([im for im in imgholder_list if im.tile_info["x"]==x_coord and im.tile_info["y"]==y_coord]))
                    counter += 1

            tile = tf.argmax(tf.reduce_mean(tile, axis=-1), axis=-1)

            with tf.device("/cpu:0"):
                tile = tile.numpy().astype("int8")
                imgholder.mask_final = tile[Start2:Start2+imageSize,Start2:Start2+imageSize]
        
        except Exception as e:
            print(f"{imgholder.filename} is fail")
            print("== error info =="); print(f"type: {str(type(e))}")
            # print(f"args: {str(e.args)}")
            print(f"error: {str(e)}")
            fail_list.append(imgholder.path)
    
    ##########
    imgholder_list = []
    print("--- Registering image info ---")
    for i in tqdm(range(len(path_list))):
        try:
            imgholder_list.append(ImageHolder(path=path_list[i]))
            
        except Exception as e:
            print(f"{path_list[i]} is fail")
            print("== error info =="); print(f"type: {str(type(e))}"); print(f"error: {str(e)}")
            fail_list.append(path_list[i])
            
    ##########
    print(f"--- Predicting images (batch size: {batch_size})---")
    for i in tqdm(range((len(imgholder_list)-1)//batch_size + 1)):
        start_i = i*batch_size; end_i = start_i + batch_size
        imgholder_batch = imgholder_list[start_i:end_i]
        img_batch = [imgholder.image for imgholder in imgholder_batch]
        img_batch = tf.stack(img_batch)
        masks = model.predict_on_batch(img_batch) 
        if softmax is True:
            masks = K.softmax(masks, axis=-1).numpy()
            
        with tf.device("/cpu:0"):
            for i in range(len(masks)):
                imgholder_batch[i].mask = masks[i]
                imgholder_batch[i].image = None
        
        #####
        for imgholder in imgholder_list:
            if imgholder.mask_final is None:
                if checkSurround(imgholder, gap) is True:
                    calcAveProb(imgholder, imageSize, overlap, numClass, gap, fail_list)
                    
        for im in imgholder_list:
            im.check_call_count()
        
    ##########
    if save_tilemask is True:
        print("--- Saving masks ---")
        with tf.device("/cpu:0"):
            for i in tqdm(range(len(imgholder_list))):
                try:
                    saveTile = imgholder_list[i].mask_final
                    saveName = "Labels_" + imgholder_list[i].filename
                    saveName = saveName[:-4]
                    saveName = saveName + ".png"
                    cv2.imwrite(f"{save_dir}/{saveName}", saveTile)    
                except:
                    pass

    if len(fail_list)>0:
        print(f"{len(fail_list)} images were fail")
        fail_list = "\n".join(fail_list)
        with open('fail_list.txt', 'a') as f:
            f.write(fail_list)
    
    return imgholder_list


def AverageProbmapSlides(
    model,
    path_list:list[str]=None,
    image_dir:str="tiles_512jpg",
    extension="jpg",
    numClass:int=14,
    imageSize:int=512,
    overlap:int=256,
    downsample:int=5,
    batch_size:int=10,
    save_tilemask:bool=False,
    save_dir:str="Pred_masks",
    model_path:str=None,
    custom_object:bool=True,
    softmax:bool=True,
    makeColorMask:bool=True,
    colormap_path="index.npy"
    ):
    """
    Runs AverageProbmap for each slide and then directly performs tiling.
    """
    
    if model is None:
        model = GetTrainedModel(model_path=model_path, custom_object=custom_object)  
        
    if path_list is None:
        imgPath = sorted(glob(f"{image_dir}/*{extension}"))
    else:
        imgPath = sorted(path_list)
    
    slidenames = sorted(GetSlideNames(path_list=imgPath,file_dir=None,prefix=None,extension=extension))
    
    global fail_list; fail_list=[]
    
    for slide in tqdm(slidenames):
        try:
            print(f"<<<<<< Predicting {slide} >>>>>>")
            paths = [s for s in imgPath if slide in s]
            print(f"{len(paths)} images exist")
            imgholder_list = AverageProbmap(
                model=model,path_list=paths,numClass=numClass, imageSize=imageSize,overlap=overlap,downsample=downsample,
                batch_size=batch_size, save_tilemask=save_tilemask, save_dir=save_dir, softmax=softmax)

            MergeTiles_from_ImageHolder(
                imgholder_list=imgholder_list, overlap_px=overlap, imageSize=imageSize, downsample=downsample,
                save_dir=save_dir, makeColorMask=makeColorMask, save_ColorMask_dir=f"{save_dir}_color", colormap_path=colormap_path)

            print(f"{slide} is done!!\n")
            K.clear_session()
        
        except Exception as e:
            print(f"{slide} is fail")
            print("== error info =="); print(f"type: {str(type(e))}")
            print(f"args: {str(e.args)}"); print(f"error: {str(e)}")
            
def SimplePrediction(path_list, 
                     save_dir="Pred_masks",
                     model_path = "Model_11class_220916.h5",
                     imageSize=256,
                     custom_object=True,
                     save_format_npy=False):
    """
    A series of functions that loads a trained model and generates predicted masks from brightfield images.
    This version simply creates mask images without considering overlapping regions.
    """
    model = GetTrainedModel(model_path=model_path, custom_object=custom_object)
    
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    for i in tqdm(range(len(image_paths))):
        path = image_paths[i]
        tmp = read_pred_image(image_path=path, imageSize=imageSize)
        tmp = model.predict_on_batch(np.expand_dims(tmp, axis=0))
        tmp = np.squeeze(tmp)
        tmp = np.argmax(tmp, axis=-1).astype("int8")
        
        saveName = os.path.basename(path)
        saveName = "Labels_" + saveName
        saveName = saveName[:-4]
        
        if save_format_npy==True:
            np.save(f"{save_dir}/{saveName}", tmp)
        else:
            saveName = saveName + ".png"
            cv2.imwrite(f"{save_dir}/{saveName}", tmp)

