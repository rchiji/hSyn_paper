## 230522 混合精度に対応させる。

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
from tensorflow.keras import mixed_precision

# initializerの引数などを事前に定義したdict。tf.keras.initializers.VarianceScalingをデフォルト値を変えて使うために予めdictで引数をまとめておく。
CONV_KERNEL_INITIALIZER = {
    "class_name": "VarianceScaling",
    "config": {"scale": 2.0,"mode": "fan_out","distribution": "truncated_normal"}
}
# --> ResNetRSではConv2Dもscale 2のVarianceScalingを使っているが、計算式上はHeNormal initializerとほぼ同じ。
# --> fan_inはinputの重みのunit数。fan_outはoutputの重みのunit数。
# --> trancated_normalは切断正規分布を使う意味。デフォルト引数。https://qiita.com/kai0706/items/bd3816342d725508fec2



## conv_blockでConv2D->BatchNormかConv2D->BatchNorm->Reluか選べるように改変
## strides=2に場合も選択できるように改変
def convolution_block(
    num_filters:int,   
    kernel_size:int=3,
    strides:int=1,
    relu:bool=True,
    dilation_rate:int=1,
    padding="same",
    use_bias=False,
):
    """""""""""""""""""""""
    convolution_block
    ・Conv2D->BatchNormかConv2D->BatchNorm->Activationの1セット
    ・relu=引数で最後に活性化関数をかけるかどうか指定
    """""""""""""""""""""""
    def apply(inputs):
        
        x = layers.Conv2D(
            filters=num_filters,
            kernel_size=kernel_size,
            strides=strides,
            dilation_rate=dilation_rate,
            padding="same",
            use_bias=use_bias,
            kernel_initializer=CONV_KERNEL_INITIALIZER, # <-- initializerをResNetRSと同様に
        )(inputs)

        x = layers.BatchNormalization()(x)

        ## relu引数がTrueの時だけ活性化関数をかける。
        if relu==True:
            x = keras.activations.relu(x)
        return x
    
    return apply


"""""""""""""""""""""""
Squeeze and Excitation block
・GlobalAveragePooling2Dでch単位の平均値を取得
・畳み込み2回 (Conv2Dの中でreluもやる。bnはなし。)
・入力データに上記の結果を行列積する。
"""""""""""""""""""""""
def SE_block(
    num_filters:int, 
    se_denomi:int = 8
):

    def apply(inputs):
        x = layers.GlobalAveragePooling2D()(inputs)
        
        # 一次元化されてしまうので、3軸のtensor形状に戻す。
        se_shape = (1, 1, x.shape[-1])
        x = layers.Reshape(se_shape)(x)

        x = layers.Conv2D(
            filters=num_filters//se_denomi,
            kernel_size=[1, 1],
            strides=[1, 1],
            kernel_initializer=CONV_KERNEL_INITIALIZER,
            padding="same",
            use_bias=True,
            activation="relu",
        )(x)

        x = layers.Conv2D(
            filters=num_filters,
            kernel_size=[1, 1],
            strides=[1, 1],
            kernel_initializer=CONV_KERNEL_INITIALIZER,
            padding="same",
            use_bias=True,
            activation="sigmoid",
        )(x)

        return layers.multiply([inputs, x])

    return apply ## <-- 返り値を関数にすることでSE(引数)(input)みたいなlayers.Conv2D(引数)(input)のような使い方ができる。



def residual_block(
    num_filters:int = 256,
    strides:int = 1,
    SE:bool = True,
    se_denomi:bool = 8,
    drop_rate:float = 0):
    """""""""""""""""""""""
    residual blockを定義
    == 畳み込みブロック ==
    ・1*1 kernelで特徴量数1/4に畳み込み
    ・3*3 kernelで特徴量数1/4に畳み込み。
    　この時strides=2にしていると画像サイズが半分にダウンサンプルできる。
    ・1*1 kernelで指定特徴量数に畳み込み
    ・SE_blockで処理
    == shortcut ブロック ==
    ・shortcutは画像サイズが半分にするときはstrides 2のAveragePoolingを通す。
    ・inputとoutputの特徴量数が違う場合は1*1 kernelの畳み込みで特徴量数を揃える。
    == add ==
    ・畳み込みブロックとshortcutブロックをaddして活性化関数を通して終了
    """""""""""""""""""""""

    def apply(inputs):
        x = convolution_block(num_filters=num_filters//4, kernel_size=1, relu=True, strides=1)(inputs) # <-- 1*1カーネルでblock_input filter数の1/4の特徴量output
        x = convolution_block(num_filters=num_filters//4, kernel_size=3, relu=True, strides=strides)(x) # <-- 3*3カーネルで同じ特徴量数のoutput
        x = convolution_block(num_filters=num_filters, kernel_size=1, relu=False, strides=1)(x) # <-- 1*1カーネルでblock_inputと同じ特徴量数のoutput
        if SE is True:
            x = SE_block(num_filters=num_filters, se_denomi = se_denomi)(x)
        
        if drop_rate > 0:
            x = layers.Dropout(rate=drop_rate, noise_shape=(None,1,1,1))(x)

        shortcut = inputs
        # ダウンサンプルがあるときはストライド2のAveragePooling -> 1*1 kernelの畳み込み
        if strides == 2:
            shortcut = layers.AveragePooling2D(pool_size=(2,2),strides=strides,padding="same")(shortcut)
            shortcut = convolution_block(num_filters=num_filters, kernel_size=1, relu=False, strides=1)(shortcut)
        
        # もし特徴量数が違っていたら1*1 kernelの畳み込み
        if shortcut.shape[-1] != x.shape[-1]:
            shortcut = convolution_block(num_filters=num_filters, kernel_size=1, relu=False, strides=1)(shortcut)
            
        out = layers.add([x, shortcut])
        out = keras.activations.relu(out)
        return out
    
    return apply # <-- 関数自体を返す。


def CustomModel_EffiNetV2(image_size, num_classes, mixed_float16=True):
    model_input = keras.Input(shape=(image_size, image_size, 3))
    
    # 混合精度の宣言
    if mixed_float16 is True:
        mixed_precision.set_global_policy('mixed_float16')
        
    EffiNetV2 = keras.applications.EfficientNetV2S(
        weights="imagenet", include_top=False, input_tensor=model_input,
        include_preprocessing = False,
    )   
    ## 512*512 inputを32*32まで落とした層を取得
    encode_16 = EffiNetV2.get_layer("top_activation").output # <-- (16,16,1280)のところ
    encode_16 = convolution_block(num_filters=512, kernel_size=1)(encode_16)
    encode_16 = residual_block(num_filters=512)(encode_16)
    encode_16 = residual_block(num_filters=512)(encode_16)      
    
    encode_8 = residual_block(num_filters=1024, strides=2)(encode_16) # <-- stride2で画像サイズを下げる。
    encode_8 = residual_block(num_filters=1024)(encode_8)
    encode_8 = residual_block(num_filters=1024)(encode_8)
    
    encode_4 = residual_block(num_filters=2048, strides=2)(encode_8)
    encode_4 = residual_block(num_filters=2048)(encode_4)
    encode_4 = residual_block(num_filters=2048)(encode_4)
    
    ## 4*4のencodeを2倍にUpsampling
    decode_8 = layers.UpSampling2D(size=(2,2), interpolation="bilinear")(encode_4)
    decode_8 = residual_block(num_filters=1024)(decode_8)
    decode_8 = residual_block(num_filters=1024)(decode_8)
    decode_8 = residual_block(num_filters=1024)(decode_8)
    
    ## 8*8のencodeを2倍にUpsampling
    decode_16 = layers.UpSampling2D(size=(2,2), interpolation="bilinear")(decode_8)
    decode_16 = residual_block(num_filters=512)(decode_16)
    decode_16 = residual_block(num_filters=512)(decode_16)
    decode_16 = residual_block(num_filters=512)(decode_16)
    
    ## 16*16のencodeを2倍にUpsampling
    decode_32 = layers.UpSampling2D(size=(2,2), interpolation="bilinear")(decode_16)
    ## 32*32のencoder-decoderをconcatinate
    encode_32 = EffiNetV2.get_layer("block6a_expand_activation").output # <-- (32,32,960)のところ
    encode_32 = convolution_block(num_filters=256, kernel_size=1)(encode_32)
    decode_32 = layers.Concatenate(axis=-1)([decode_32, encode_32])    
    decode_32 = residual_block(num_filters=256)(decode_32) 
    decode_32 = residual_block(num_filters=256)(decode_32) 
    decode_32 = residual_block(num_filters=256)(decode_32) 
    
    ## 2倍にupsampling
    decode_64 = layers.UpSampling2D(size=(2,2), interpolation="bilinear")(decode_32)
    # encoderの64*64層を取得
    encode_64 = EffiNetV2.get_layer("block4a_expand_activation").output # <-- (64,64,256)のところ
    encode_64 = convolution_block(num_filters=128, kernel_size=1)(encode_64)
    decode_64 = layers.Concatenate(axis=-1)([decode_64, encode_64]) 
    decode_64 = residual_block(num_filters=128)(decode_64)
    decode_64 = residual_block(num_filters=128)(decode_64)
    
    ## 2倍にupsampling
    decode_128 = layers.UpSampling2D(size=(2,2), interpolation="bilinear")(decode_64)
    # encoderの128*128層を取得
    encode_128 = EffiNetV2.get_layer("block2d_add").output # <-- (128,128,48)のところ
    encode_128 = convolution_block(num_filters=64, kernel_size=1)(encode_128)
    decode_128 = layers.Concatenate(axis=-1)([decode_128, encode_128]) 
    decode_128 = residual_block(num_filters=64)(decode_128)
    decode_128 = residual_block(num_filters=64)(decode_128)   
    
    ## 4倍にUpsamplingして元の512*512に
    x = layers.UpSampling2D(size=(4,4), interpolation="bilinear")(decode_128)
    
    model_output = layers.Conv2D(num_classes, kernel_size=(1, 1), padding="same", dtype="float32")(x)
    return keras.Model(inputs=model_input, outputs=model_output)


def CustomModel_ResNetRS101(image_size, num_classes, mixed_float16=True):
    model_input = keras.Input(shape=(image_size, image_size, 3))
    
    # 混合精度の宣言
    if mixed_float16 is True:
        mixed_precision.set_global_policy('mixed_float16')
    
    resnet_rs101 = keras.applications.resnet_rs.ResNetRS101(
        weights="imagenet", include_top=False, input_tensor=model_input,
        include_preprocessing = False,
    )


    ## 512*512 inputを32*32まで落とした層を取得
    encode_32 = resnet_rs101.get_layer("BlockGroup4__block_22__output_act").output # <-- (32,32,1024)のところ
    
    encode_16 = residual_block(num_filters=512, strides=2, drop_rate=0.24)(encode_32) # <-- stride2で画像サイズを下げる。
    encode_16 = residual_block(num_filters=512, drop_rate=0.24)(encode_16)
    encode_16 = residual_block(num_filters=512, drop_rate=0.24)(encode_16)   
    
    encode_8 = residual_block(num_filters=1024, strides=2, drop_rate=0.28)(encode_16) # <-- stride2で画像サイズを下げる。
    encode_8 = residual_block(num_filters=1024, drop_rate=0.28)(encode_8)
    encode_8 = residual_block(num_filters=1024, drop_rate=0.28)(encode_8)
    
    encode_4 = residual_block(num_filters=2048, strides=2, drop_rate=0.32)(encode_8)
    encode_4 = residual_block(num_filters=2048, drop_rate=0.32)(encode_4)
    encode_4 = residual_block(num_filters=2048, drop_rate=0.32)(encode_4)
    
    ## 4*4のencodeを2倍にUpsampling
    decode_8 = layers.UpSampling2D(size=(2,2), interpolation="bilinear")(encode_4)
    for i in range(3):
        decode_8 = residual_block(num_filters=1024, SE=False)(decode_8)
    
    ## 8*8のencodeを2倍にUpsampling
    decode_16 = layers.UpSampling2D(size=(2,2), interpolation="bilinear")(decode_8)
    for i in range(3):
        decode_16 = residual_block(num_filters=512, SE=False)(decode_16)
    
    ## 16*16のencodeを2倍にUpsampling
    decode_32 = layers.UpSampling2D(size=(2,2), interpolation="bilinear")(decode_16)
    ## 32*32のencoder-decoderをconcatinate
    encode_32 = convolution_block(encode_32, num_filters=48, kernel_size=1)
    decode_32 = layers.Concatenate(axis=-1)([decode_32, encode_32])    
    for i in range(3):
        decode_32 = residual_block(num_filters=256, SE=False)(decode_32)
    
    ## 2倍にupsamplingして64*64にしてdsppを通す
    decode_64 = layers.UpSampling2D(size=(2,2), interpolation="bilinear")(decode_32)    
    decode_64 = residual_block(num_filters=128, SE=False)(decode_64)
    decode_64 = residual_block(num_filters=64, SE=False)(decode_64)
    
    ## 8倍にUpsamplingして元の512*512に
    x = layers.UpSampling2D(size=(8,8), interpolation="bilinear")(decode_64)
    
    model_output = layers.Conv2D(num_classes, kernel_size=(1, 1), padding="same", dtype="float32")(x)
    return keras.Model(inputs=model_input, outputs=model_output)