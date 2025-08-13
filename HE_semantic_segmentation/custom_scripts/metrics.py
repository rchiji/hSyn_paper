import os
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'
import numpy as np
import tensorflow as tf

from losses import *

def precision(y_true, y_pred):
    smooth = 1e-7
    # 正解ラベル、予測ラベルの前処理
    y_true,y_pred = preprocess(y_true,y_pred)
    
    # score calculation
    tp = y_true * y_pred
    fp = y_pred - tp

    tp = tf.reduce_sum(tp,axis=[0,1,2])
    fp = tf.reduce_sum(fp, axis=[0,1,2])
    
    score = (tp + smooth) / (tp + fp + smooth)
    
    return tf.reduce_mean(score)

def recall(y_true, y_pred):
    smooth = 1e-7
    # 正解ラベル、予測ラベルの前処理
    y_true,y_pred = preprocess(y_true,y_pred)
    
    # score calculation
    tp = y_true * y_pred
    fn = y_true - tp

    tp = tf.reduce_sum(tp,axis=[0,1,2])
    fn = tf.reduce_sum(fn,axis=[0,1,2])
    
    score = (tp + smooth) / (tp + fn + smooth)
    
    return tf.reduce_mean(score)


def dice(y_true, y_pred):
    smooth=1e-7
    # 正解ラベル、予測ラベルの前処理
    y_true,y_pred = preprocess(y_true,y_pred)

    # [Batch*Y*X, num_classes]の形状に
    y_true = tf.reshape(y_true, [-1,tf.shape(y_true)[-1]])
    y_pred = tf.reshape(y_pred, [-1,tf.shape(y_pred)[-1]])

    # diceの分子、分母を計算 [Batch*Y*X, num_classes] -> [num_classes]
    intersect = tf.reduce_sum(y_true*y_pred, axis=0)
    denom = tf.reduce_sum(y_true+y_pred, axis=0)

    dice_values = (2.*intersect+smooth)/(denom+smooth)
    return tf.reduce_mean(dice_values)


def jaccard(y_true, y_pred):
    smooth = 1e-7   
    # 正解ラベル、予測ラベルの前処理
    y_true,y_pred = preprocess(y_true,y_pred)

    # [Batch*Y*X, num_classes]の形状に
    y_true = tf.reshape(y_true, [-1,tf.shape(y_true)[-1]])  
    y_pred = tf.reshape(y_pred, [-1,tf.shape(y_pred)[-1]])
    # jaccardの分子、分母を計算
    intersect = tf.reduce_sum(y_true*y_pred, axis=0)
    union = tf.reduce_sum(y_true+y_pred-(y_true*y_pred), axis=0) 
    jaccard = (intersect + smooth) / (union + smooth)

    return tf.reduce_mean(jaccard)


class Dice(tf.keras.metrics.Metric):
    def __init__(self, 
                 class_average=True,
                 smooth=1e-7,
                 name="dice",
                 **kwargs):
        super().__init__(name=name, **kwargs)
        self.class_average=class_average
        self.smooth=smooth
        self.mean = tf.metrics.Mean()

    # stepの度に呼ばれる関数。評価値を計算する
    def update_state(self, y_true, y_pred, sample_weight=None):
        # 分類class数を自動判定
        num_classes = tf.cast(tf.shape(y_pred)[-1],tf.int32)

        # 正解ラベル、予測ラベルの前処理
        y_true,y_pred = preprocess(y_true,y_pred)
    
        # 平均する軸以外の次元をまとめる
        shape = [-1,num_classes] if self.class_average else [-1]

        y_true_f = tf.reshape(y_true, shape)
        y_pred_f = tf.reshape(y_pred, shape)
    
        # diceの分子、分母を計算
        intersect = tf.reduce_sum(y_true_f * y_pred_f,axis=0)
        denom = tf.reduce_sum(y_true_f + y_pred_f,axis=0)
        # diceの値
        values = (2. * intersect + self.smooth) / (denom + self.smooth)
        
        if sample_weight is not None:
            sample_weight = tf.cast(sample_weight, "float32")
            values = tf.multiply(values, sample_weight)
        
        # 結果を登録
        self.mean.update_state(values)
        
    # stepの度に呼ばれる関数。結果を表示する。
    def result(self):
        return self.mean.result()

    # epochの度に呼ばれる関数。値をリセットする。
    def reset_state(self):
        self.mean.reset_state()


class Dice_perclass(tf.keras.metrics.Metric):
    def __init__(self, 
                 smooth=1e-7,
                 name="dice_perclass",
                 **kwargs):
        super().__init__(name=name, **kwargs)
        self.smooth=smooth
        self.value = None
        self.num_classes = None

    def update_state(self, y_true, y_pred, sample_weight=None):    
        
        if self.value is None:
            # 分類class数を自動判定
            self.num_classes = y_pred.shape.as_list()[-1]
            self.value = self.add_weight(name=self.name,
                                         shape=(self.num_classes),
                                         initializer="zeros")

        # 正解ラベル、予測ラベルの前処理
        y_true,y_pred = preprocess(y_true,y_pred)
    
        # 平均する軸以外の次元をまとめる
        shape = [-1, self.num_classes]
        y_true_f = tf.reshape(y_true, shape)
        y_pred_f = tf.reshape(y_pred, shape)
    
        # diceの分子、分母を計算
        intersect = tf.reduce_sum(y_true_f * y_pred_f, axis=0)
        denom = tf.reduce_sum(y_true_f + y_pred_f, axis=0)
        # diceの値
        values = (2. * intersect + self.smooth) / (denom + self.smooth)

        if sample_weight is not None:
            sample_weight = tf.cast(sample_weight, "float32")
            values = tf.multiply(values, sample_weight)
        
        # 結果を登録
        self.value.assign(values)
        
    def result(self):
        return self.value

    def reset_state(self):
        self.value.assign(tf.zeros_like(self.value))