## 221031 update

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
    # one_hot化したtrue_labelとの積なのでtrue_labelの予測確率値のみが残る。
    true_posi = K.sum(y_true_f * y_pred_f, axis=-1) # --> TP
    # 正解posi数, 予測posi数をすべて足し合わせると2TP+FP+FNになるので、TPを一つだけ引く
    denom = K.sum(y_true_f + y_pred_f, axis=-1) - true_posi # --> TP+FP+FN
    return K.mean((true_posi / (denom + smooth)))

def IoU_loss(y_true, y_pred):
    return 1 - IoU(y_true, y_pred,smooth=1e-7)

def Dice_SCCE_loss(y_true, y_pred, weight:int=2,from_logits:bool=True):
    
    custom_loss = Dice_loss(y_true, y_pred)
    
    # loss関数をオブジェクト化して使用
    scce = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=from_logits)
    scce = scce(y_true, y_pred)
     
    return custom_loss*weight + scce

def IoU_SCCE_loss(y_true, y_pred, weight:int=2,from_logits:bool=True):
    
    custom_loss = IoU_loss(y_true, y_pred)
    
    # loss関数をオブジェクト化して使用
    scce = tf.keras.losses.SparseCategoricalCrossentropy(from_logits=from_logits)
    scce = scce(y_true, y_pred)
     
    return custom_loss*weight + scce

"""
ネットに落ちてたUpadatedMeanIoU
https://stackoverflow.com/questions/61824470/dimensions-mismatch-error-when-using-tf-metrics-meaniou-with-sparsecategorical
https://github.com/tensorflow/tensorflow/issues/32875
tf.keras.metrics.MeanIoUの親クラスのupdate_stateで混合行列作ったりするのだが、マルチクラス用にargmaxを入れてからupdate_state関数が進むように改変したもの。
"""
class UpdatedMeanIoU(tf.keras.metrics.MeanIoU):
    def __init__(self,
                 y_true=None,
                 y_pred=None,
                 num_classes=None,
                 name="MeanIoU",
                 dtype=None):
        super().__init__(num_classes = num_classes,name=name, dtype=dtype)
    
    # 親クラスに既にあるupdate_state関数をオーバーライド。argmaxするように改変。
    def update_state(self, y_true, y_pred, sample_weight=None):
        y_pred = tf.math.argmax(y_pred, axis=-1)
        
        return super().update_state(y_true, y_pred, sample_weight)

"""
tf.keras.metrics.MeanIoUのsource codeを参照して書いたdef
update_state関数、result関数と同様の構造を取っているのだが、なぜかMeanIoUクラスを引き継いだ上記のモデルの再現ができない。
"""
def MeanIoU2(y_true, y_pred, sample_weight=None):
    num_classes = y_pred.shape[-1]
    # 予測値はクラスの数だけあるので、argmaxする。
    y_pred = tf.argmax(y_pred, axis=-1)
    # データ型をデフォルトにしておく
    y_true = tf.cast(y_true, K.floatx())
    y_pred = tf.cast(y_pred, K.floatx())

    # 軸数が1より大きければ1次元化する。(batch * widht * heightの要素数)
    if y_pred.shape.ndims > 1:
        y_pred = tf.reshape(y_pred, [-1])
    if y_true.shape.ndims > 1:
        y_true = tf.reshape(y_true, [-1])

    if sample_weight is not None:
        sample_weight = tf.cast(sample_weight, K.floatx())
        if sample_weight.shape.ndims > 1:
            sample_weight = tf.reshape(sample_weight, [-1])

    # 混合行列作成 0,1
    cm = tf.math.confusion_matrix(labels=y_true,predictions=y_pred,num_classes=num_classes,
                                  weights=sample_weight, dtype=K.floatx())
    # 行方向の合計
    sum_over_row = tf.cast(tf.reduce_sum(cm, axis=0),dtype=K.floatx())
    # 列方向の合計
    sum_over_col = tf.cast(tf.reduce_sum(cm, axis=1),dtype=K.floatx())
    # tensorの対角成分(TP)の合計
    true_positives = tf.cast(tf.linalg.tensor_diag_part(cm),dtype=K.floatx())
    
    # 分母 (TP+FP+FN)の作成
    denominator = sum_over_row + sum_over_col - true_positives # 行方向、列方向で対角成分が2回重複して数えられているので1回分TPを引く

    # Only keep the target classes <-- このblockは無くていいはずやけど。特定のclassだけを評価させることもできそう。
    target_class_ids=list(range(num_classes))
    true_positives = tf.gather(true_positives, target_class_ids)
    denominator = tf.gather(denominator, target_class_ids)

    # If the denominator is 0, we need to ignore the class.
    num_valid_entries = tf.reduce_sum(tf.cast(tf.not_equal(denominator, 0), dtype=K.floatx()))

    iou = tf.math.divide_no_nan(true_positives, denominator)

    return tf.math.divide_no_nan(tf.reduce_sum(iou), num_valid_entries)