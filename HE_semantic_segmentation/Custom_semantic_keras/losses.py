import os
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'
import re
import random
import numpy as np
from glob import glob
from scipy.io import loadmat
import matplotlib.pyplot as plt
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import keras.backend as K

def Dice(y_true, y_pred,smooth=1e-7):
    num_classes = y_pred.shape[-1]
    y_pred = K.softmax(y_pred)
    y_true_f = K.flatten(K.one_hot(K.cast(y_true, 'int32'), num_classes=num_classes))
    y_pred_f = K.flatten(y_pred)
    intersect = K.sum(y_true_f * y_pred_f, axis=-1)
    denom = K.sum(y_true_f + y_pred_f, axis=-1)
    return K.mean((2. * intersect / (denom + smooth)))

def Dice_loss(y_true, y_pred):
    return 1 - Dice(y_true, y_pred,smooth=1e-7)

def IoU(y_true, y_pred,smooth=1e-7):
    num_classes = y_pred.shape[-1]
    y_pred = K.softmax(y_pred)
    y_true_f = K.flatten(K.one_hot(K.cast(y_true, 'int32'), num_classes=num_classes))
    y_pred_f = K.flatten(y_pred)
    true_posi = K.sum(y_true_f * y_pred_f, axis=-1) # --> TP
    denom = K.sum(y_true_f + y_pred_f, axis=-1) - true_posi # --> TP+FP+FN
    return K.mean((true_posi / (denom + smooth)))

def IoU_loss(y_true, y_pred):
    return 1 - IoU(y_true, y_pred,smooth=1e-7)

def Dice_SCCE_loss(y_true, y_pred, weight:int=2,from_logits:bool=True):
    
    custom_loss = Dice_loss(y_true, y_pred)
    
    scce = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=from_logits)
    scce = scce(y_true, y_pred)
     
    return custom_loss*weight + scce

def IoU_SCCE_loss(y_true, y_pred, weight:int=2,from_logits:bool=True):
    
    custom_loss = IoU_loss(y_true, y_pred)
    
    scce = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=from_logits)
    scce = scce(y_true, y_pred)
     
    return custom_loss*weight + scce


class UpdatedMeanIoU(tf.keras.metrics.MeanIoU):
    def __init__(self,
                 y_true=None,
                 y_pred=None,
                 num_classes=None,
                 name="MeanIoU",
                 dtype=None):
        super().__init__(num_classes = num_classes,name=name, dtype=dtype)
    
    def update_state(self, y_true, y_pred, sample_weight=None):
        y_pred = tf.math.argmax(y_pred, axis=-1)
        
        return super().update_state(y_true, y_pred, sample_weight)


def MeanIoU2(y_true, y_pred, sample_weight=None):
    num_classes = y_pred.shape[-1]
    y_pred = tf.argmax(y_pred, axis=-1)
    y_true = tf.cast(y_true, K.floatx())
    y_pred = tf.cast(y_pred, K.floatx())

    if y_pred.shape.ndims > 1:
        y_pred = tf.reshape(y_pred, [-1])
    if y_true.shape.ndims > 1:
        y_true = tf.reshape(y_true, [-1])

    if sample_weight is not None:
        sample_weight = tf.cast(sample_weight, K.floatx())
        if sample_weight.shape.ndims > 1:
            sample_weight = tf.reshape(sample_weight, [-1])

    cm = tf.math.confusion_matrix(labels=y_true,predictions=y_pred,num_classes=num_classes,
                                  weights=sample_weight, dtype=K.floatx())
    sum_over_row = tf.cast(tf.reduce_sum(cm, axis=0),dtype=K.floatx())
    sum_over_col = tf.cast(tf.reduce_sum(cm, axis=1),dtype=K.floatx())
    true_positives = tf.cast(tf.linalg.tensor_diag_part(cm),dtype=K.floatx())
    
    denominator = sum_over_row + sum_over_col - true_positives

    target_class_ids=list(range(num_classes))
    true_positives = tf.gather(true_positives, target_class_ids)
    denominator = tf.gather(denominator, target_class_ids)

    num_valid_entries = tf.reduce_sum(tf.cast(tf.not_equal(denominator, 0), dtype=K.floatx()))

    iou = tf.math.divide_no_nan(true_positives, denominator)

    return tf.math.divide_no_nan(tf.reduce_sum(iou), num_valid_entries)