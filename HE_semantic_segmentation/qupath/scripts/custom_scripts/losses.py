import os
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'
import numpy as np
import tensorflow as tf
from tqdm.auto import tqdm

from datasets import prepare_image_paths


def preprocess(y_true,y_pred):
    # 分類class数を自動判定
    num_classes = tf.shape(y_pred)[-1]

    # SoftmaxされていなければSoftmax
    if not tf.round(tf.reduce_mean(tf.reduce_sum(y_pred[0],axis=-1))*10)/10 == 1.0:
        y_pred = tf.nn.softmax(y_pred, axis=-1)
    
    # 正解ラベルをone hot encoding
    if len(tf.shape(y_true)) == 3:
        y_true = tf.one_hot(tf.cast(y_true, dtype=tf.uint8), depth=num_classes)

    return y_true,y_pred


def dice_loss(y_true, y_pred):
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
    
    return 1.0 - tf.reduce_mean(dice_values)


def jaccard_loss(y_true, y_pred):
    smooth = 1e-7
    
    # 正解ラベル、予測ラベルの前処理
    y_true,y_pred = preprocess(y_true,y_pred)

    # [Batch*Y*X, num_classes]の形状に
    y_true = tf.reshape(y_true, [-1,tf.shape(y_true)[-1]])
    y_pred = tf.reshape(y_pred, [-1,tf.shape(y_pred)[-1]])
    
    # jaccardの分子、分母を計算 [Batch*Y*X, num_classes] -> [num_classes]
    intersect = tf.reduce_sum(y_true * y_pred, axis=0)
    union = tf.reduce_sum(y_true + y_pred - (y_true * y_pred), axis=0) 
    jaccard = (intersect + smooth) / (union + smooth)

    return 1.0 - tf.reduce_mean(jaccard)
    
    
class DiceLoss(tf.keras.losses.Loss):
    def __init__(self, 
                 class_average=True,
                 smooth=1e-7,
                 name="dice_loss",
                 class_weight=None
                ):
        super().__init__(name=name)
        self.class_average=class_average
        self.smooth=smooth
        self.class_weight = class_weight

    def call(self, y_true, y_pred):
        # 正解ラベル、予測ラベルの前処理
        y_true,y_pred = preprocess(y_true,y_pred)
    
        # 形状変化。class_average=Trueなら[Batch*Y*X, num_classes], Falseなら[Batch*Y*X*num_classes]
        shape = [-1,tf.shape(y_true)[-1]] if self.class_average is True else [-1]
        y_true_f = tf.reshape(y_true, shape)
        y_pred_f = tf.reshape(y_pred, shape)
    
        # diceの分子、分母を計算
        intersect = tf.reduce_sum(y_true_f * y_pred_f, axis=0)
        denom = tf.reduce_sum(y_true_f + y_pred_f, axis=0)
    
        dice_values = (2. * intersect + self.smooth) / (denom + self.smooth)
        dice_loss = 1.0 - dice_values

        if self.class_average and self.class_weight is not None:
            dice_loss = dice_loss * self.class_weight
            
        return tf.reduce_mean(dice_loss)


class JaccardLoss(tf.keras.losses.Loss):

    def __init__(self,
                 class_average=True, 
                 smooth=1e-7, 
                 name="jaccard_loss"):
        super().__init__(name=name)
        self.class_average = class_average 
        self.smooth = smooth

    def call(self, y_true, y_pred):
        # 正解ラベル、予測ラベルの前処理
        y_true,y_pred = preprocess(y_true,y_pred)
    
        # 形状変化。class_average=Trueなら[Batch*Y*X, num_classes], Falseなら[Batch*Y*X*num_classes]
        shape = [-1,tf.shape(y_true)[-1]] if self.class_average is True else [-1]
        y_true_f = tf.reshape(y_true, shape)
        y_pred_f = tf.reshape(y_pred, shape)

        # jaccardの分子、分母を計算
        intersect = tf.reduce_sum(y_true_f * y_pred_f, axis=0)
        union = tf.reduce_sum(y_true_f + y_pred_f - (y_true_f * y_pred_f), axis=0)
        jaccard = (intersect + self.smooth) / (union + self.smooth)  

        return 1.0 - tf.reduce_mean(jaccard)


def CalcLossWeight(
    label_dir:str,
    num_classes:int,
    log_base=None
):
    # ラベル画像のファイルパスリスト作成
    label_paths = prepare_image_paths(label_dir)

    # ラベル画像を読み込んで各ラベル輝度のピクセル数を集計
    index_counts = []
    for path in tqdm(label_paths):
        # ラベル画像読み込み
        label = tf.io.read_file(path)
        label = tf.image.decode_png(label, channels=1)
        # ラベル輝度ごとに集計
        index_count = []
        for i in range(num_classes):
            index_count.append(np.count_nonzero(label == i))
        index_counts.append(index_count)

    class_weight = np.sum(np.array(index_counts))/np.sum(np.array(index_counts),axis=0)
    
    if log_base is not None:
        import math
        class_weight = np.array([math.log(cw,log_base) for cw in class_weight])
    
    return class_weight


class SSCEDiceLoss(tf.keras.losses.Loss):
    def __init__(self, 
                 loss_weight=[1,1],
                 class_average=True, 
                 from_logits=False,
                 class_weight=None,
                 smooth=1e-7,
                 name="ssce_dice_loss"):
        super().__init__(name=name)
        self.dice_loss=DiceLoss(
            class_average=class_average,
            class_weight=class_weight,
            smooth=smooth)
        if class_weight is not None:
            self.scce = WeightedSCCE(
                class_weight=class_weight,
                from_logits=from_logits
            )
        else:         
            self.scce = tf.keras.losses.SparseCategoricalCrossentropy(
                from_logits=from_logits
            )
        self.loss_weight=loss_weight

    def call(self, y_true, y_pred):
        scce_loss = self.scce(y_true, y_pred)
        dice_loss = self.dice_loss(y_true, y_pred)  
        return scce_loss*self.loss_weight[0] + dice_loss*self.loss_weight[1]
        

class WeightedSCCE(tf.keras.losses.Loss):
    def __init__(self, class_weight, from_logits=False, name='weighted_scce'):
        super().__init__(name=name)
        if class_weight is None or all(v == 1. for v in class_weight):
            self.class_weight = None
        else:
            self.class_weight = tf.convert_to_tensor(class_weight,
                dtype=tf.float32)
        self.unreduced_scce = tf.keras.losses.SparseCategoricalCrossentropy(
            from_logits=from_logits, name=name,
            reduction=tf.keras.losses.Reduction.NONE)

    def call(self, y_true, y_pred):
        loss = self.unreduced_scce(y_true, y_pred)
        if self.class_weight is not None:
            weight_mask = tf.gather(self.class_weight, indices=tf.cast(y_true,tf.int32))
            loss = tf.math.multiply(loss, weight_mask)
        return loss



        
