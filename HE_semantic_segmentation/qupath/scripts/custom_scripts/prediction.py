import os
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'
import numpy as np
from glob import glob
from tqdm.auto import tqdm
import tensorflow as tf
import matplotlib.pyplot as plt
from tiling import *
import cv2
import re

from datasets import *

    
def predictDataset_labelholder(
    model,
    dataset,
):
    labelholders = []
    
    for images, paths in tqdm(dataset):
        labels = model.predict_on_batch(images)
        
        for path, label in zip(paths, labels):
            labelholder = LabelHolder(path.numpy().decode())
            labelholder.label = label
            labelholders.append(labelholder)

    return labelholders

def aveprob_labelholder(
    labelholders,
    tile_size,
    overlap,
    gap,
    num_classes
):
    for labelholder in labelholders:
        if labelholder.label_final is None:
            if checkSurround(labelholder, labelholders, gap) is True:
                calcAveProb(labelholder, labelholders, tile_size,
                            overlap, gap, num_classes)
                
    for lh in labelholders:
        lh.check_call_count()

def registerLabel(
    path,
    label,
    labelholders
):
    if tf.is_tensor(path):
        path = path.numpy().decode()
    labelholder = [ lh for lh in labelholders if lh.path == path][0]
    labelholder.label = label

    
def orderLabelHolders(
    labelholders
):
    x_list = sorted(set([ lh.tile_info["x"] for lh in labelholders]))
    y_list = sorted(set([ lh.tile_info["y"] for lh in labelholders]))

    x_sets = [x_list[i:i+3] for i in range(0, len(x_list), 3)]
    y_sets = [y_list[i:i+3] for i in range(0, len(y_list), 3)]

    ordered_labelholders = []

    progressbar = tqdm(total = len(labelholders), ncols=0)
    for x_set in x_sets:
        for y_set in y_sets:
            for x in x_set:
                for y in y_set:
                    labelholder = next(iter([ lh for lh in labelholders if lh.tile_info["x"]==x and lh.tile_info["y"]==y]), None)
                    if labelholder:
                        ordered_labelholders.append(labelholder)
                    progressbar.update()
    return ordered_labelholders

def orderImagePath(
    image_paths
):
    coord_dict = { p: getTileCoord(p) for p in image_paths}
    x_list = sorted(set([ tile_info["x"] for tile_info in coord_dict.values() ]))
    y_list = sorted(set([ tile_info["y"] for tile_info in coord_dict.values() ]))

    x_sets = [x_list[i:i+3] for i in range(0, len(x_list), 3)]
    y_sets = [y_list[i:i+3] for i in range(0, len(y_list), 3)]

    ordered_path = []

    progressbar = tqdm(total = len(coord_dict))
    for x_set in x_sets:
        for y_set in y_sets:
            for x in x_set:
                for y in y_set:
                    path = next(iter([ p for p,tile_info in coord_dict.items() if tile_info["x"]==x and tile_info["y"]==y]), None)
                    if path:
                        ordered_path.append(path)
                    progressbar.update()
    return ordered_path
