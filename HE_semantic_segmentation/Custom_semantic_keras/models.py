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

CONV_KERNEL_INITIALIZER = {
    "class_name": "VarianceScaling",
    "config": {"scale": 2.0,"mode": "fan_out","distribution": "truncated_normal"}
}


def convolution_block(
    num_filters:int,   
    kernel_size:int=3,
    strides:int=1,
    relu:bool=True,
    dilation_rate:int=1,
    padding="same",
    use_bias=False,
):
    """
    convolution_block
    - One block consisting of Conv2D -> BatchNorm or Conv2D -> BatchNorm -> Activation
    - relu: specifies whether to apply an activation function at the end
    """
    def apply(inputs):
        
        x = layers.Conv2D(
            filters=num_filters,
            kernel_size=kernel_size,
            strides=strides,
            dilation_rate=dilation_rate,
            padding="same",
            use_bias=use_bias,
            kernel_initializer=CONV_KERNEL_INITIALIZER,
        )(inputs)

        x = layers.BatchNormalization()(x)

        if relu==True:
            x = keras.activations.relu(x)
        return x
    
    return apply


"""
Squeeze and Excitation block
- Computes channel-wise averages using GlobalAveragePooling2D
- Applies two convolution layers (activation included in Conv2D, no batch normalization)
- Multiplies the result with the input feature map
"""
def SE_block(
    num_filters:int, 
    se_denomi:int = 8
):

    def apply(inputs):
        x = layers.GlobalAveragePooling2D()(inputs)
        
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

    return apply 


def residual_block(
    num_filters:int = 256,
    strides:int = 1,
    SE:bool = True,
    se_denomi:bool = 8,
    drop_rate:float = 0):
    """
    Defines a residual block
    == Convolution block ==
    - Applies a 1×1 convolution to reduce the number of channels to 1/4
    - Applies a 3×3 convolution with 1/4 channels
      (if strides=2, the spatial resolution is downsampled by half)
    - Applies a 1×1 convolution to restore the specified number of channels
    - Processes the output with an SE block
    == Shortcut block ==
    - If downsampling is required, applies AveragePooling with stride 2
    - If the number of channels differs between input and output,
      applies a 1×1 convolution to match the channel dimensions
    == Add ==
    - Adds the outputs of the convolution block and shortcut block,
      then applies an activation function
    """

    def apply(inputs):
        x = convolution_block(num_filters=num_filters//4, kernel_size=1, relu=True, strides=1)(inputs)
        x = convolution_block(num_filters=num_filters//4, kernel_size=3, relu=True, strides=strides)(x)
        x = convolution_block(num_filters=num_filters, kernel_size=1, relu=False, strides=1)(x)
        if SE is True:
            x = SE_block(num_filters=num_filters, se_denomi = se_denomi)(x)
        
        if drop_rate > 0:
            x = layers.Dropout(rate=drop_rate, noise_shape=(None,1,1,1))(x)

        shortcut = inputs
        if strides == 2:
            shortcut = layers.AveragePooling2D(pool_size=(2,2),strides=strides,padding="same")(shortcut)
            shortcut = convolution_block(num_filters=num_filters, kernel_size=1, relu=False, strides=1)(shortcut)
        
        if shortcut.shape[-1] != x.shape[-1]:
            shortcut = convolution_block(num_filters=num_filters, kernel_size=1, relu=False, strides=1)(shortcut)
            
        out = layers.add([x, shortcut])
        out = keras.activations.relu(out)
        return out
    
    return apply


def CustomModel_EffiNetV2(image_size, num_classes, mixed_float16=True):
    model_input = keras.Input(shape=(image_size, image_size, 3))
    
    if mixed_float16 is True:
        mixed_precision.set_global_policy('mixed_float16')
        
    EffiNetV2 = keras.applications.EfficientNetV2S(
        weights="imagenet", include_top=False, input_tensor=model_input,
        include_preprocessing = False,
    )   
    encode_16 = EffiNetV2.get_layer("top_activation").output
    encode_16 = convolution_block(num_filters=512, kernel_size=1)(encode_16)
    encode_16 = residual_block(num_filters=512)(encode_16)
    encode_16 = residual_block(num_filters=512)(encode_16)      
    
    encode_8 = residual_block(num_filters=1024, strides=2)(encode_16)
    encode_8 = residual_block(num_filters=1024)(encode_8)
    encode_8 = residual_block(num_filters=1024)(encode_8)
    
    encode_4 = residual_block(num_filters=2048, strides=2)(encode_8)
    encode_4 = residual_block(num_filters=2048)(encode_4)
    encode_4 = residual_block(num_filters=2048)(encode_4)
    
    decode_8 = layers.UpSampling2D(size=(2,2), interpolation="bilinear")(encode_4)
    decode_8 = residual_block(num_filters=1024)(decode_8)
    decode_8 = residual_block(num_filters=1024)(decode_8)
    decode_8 = residual_block(num_filters=1024)(decode_8)
    
    decode_16 = layers.UpSampling2D(size=(2,2), interpolation="bilinear")(decode_8)
    decode_16 = residual_block(num_filters=512)(decode_16)
    decode_16 = residual_block(num_filters=512)(decode_16)
    decode_16 = residual_block(num_filters=512)(decode_16)
    
    decode_32 = layers.UpSampling2D(size=(2,2), interpolation="bilinear")(decode_16)
    encode_32 = EffiNetV2.get_layer("block6a_expand_activation").output
    encode_32 = convolution_block(num_filters=256, kernel_size=1)(encode_32)
    decode_32 = layers.Concatenate(axis=-1)([decode_32, encode_32])    
    decode_32 = residual_block(num_filters=256)(decode_32) 
    decode_32 = residual_block(num_filters=256)(decode_32) 
    decode_32 = residual_block(num_filters=256)(decode_32) 
    
    decode_64 = layers.UpSampling2D(size=(2,2), interpolation="bilinear")(decode_32)
    encode_64 = EffiNetV2.get_layer("block4a_expand_activation").output
    encode_64 = convolution_block(num_filters=128, kernel_size=1)(encode_64)
    decode_64 = layers.Concatenate(axis=-1)([decode_64, encode_64]) 
    decode_64 = residual_block(num_filters=128)(decode_64)
    decode_64 = residual_block(num_filters=128)(decode_64)
    
    decode_128 = layers.UpSampling2D(size=(2,2), interpolation="bilinear")(decode_64)
    encode_128 = EffiNetV2.get_layer("block2d_add").output
    encode_128 = convolution_block(num_filters=64, kernel_size=1)(encode_128)
    decode_128 = layers.Concatenate(axis=-1)([decode_128, encode_128]) 
    decode_128 = residual_block(num_filters=64)(decode_128)
    decode_128 = residual_block(num_filters=64)(decode_128)   
    
    x = layers.UpSampling2D(size=(4,4), interpolation="bilinear")(decode_128)
    
    model_output = layers.Conv2D(num_classes, kernel_size=(1, 1), padding="same", dtype="float32")(x)
    return keras.Model(inputs=model_input, outputs=model_output)


def CustomModel_ResNetRS101(image_size, num_classes, mixed_float16=True):
    model_input = keras.Input(shape=(image_size, image_size, 3))
    
    if mixed_float16 is True:
        mixed_precision.set_global_policy('mixed_float16')
    
    resnet_rs101 = keras.applications.resnet_rs.ResNetRS101(
        weights="imagenet", include_top=False, input_tensor=model_input,
        include_preprocessing = False,
    )


    encode_32 = resnet_rs101.get_layer("BlockGroup4__block_22__output_act").output
    
    encode_16 = residual_block(num_filters=512, strides=2, drop_rate=0.24)(encode_32)
    encode_16 = residual_block(num_filters=512, drop_rate=0.24)(encode_16)
    encode_16 = residual_block(num_filters=512, drop_rate=0.24)(encode_16)   
    
    encode_8 = residual_block(num_filters=1024, strides=2, drop_rate=0.28)(encode_16)
    encode_8 = residual_block(num_filters=1024, drop_rate=0.28)(encode_8)
    encode_8 = residual_block(num_filters=1024, drop_rate=0.28)(encode_8)
    
    encode_4 = residual_block(num_filters=2048, strides=2, drop_rate=0.32)(encode_8)
    encode_4 = residual_block(num_filters=2048, drop_rate=0.32)(encode_4)
    encode_4 = residual_block(num_filters=2048, drop_rate=0.32)(encode_4)
    
    decode_8 = layers.UpSampling2D(size=(2,2), interpolation="bilinear")(encode_4)
    for i in range(3):
        decode_8 = residual_block(num_filters=1024, SE=False)(decode_8)
    
    decode_16 = layers.UpSampling2D(size=(2,2), interpolation="bilinear")(decode_8)
    for i in range(3):
        decode_16 = residual_block(num_filters=512, SE=False)(decode_16)
    
    decode_32 = layers.UpSampling2D(size=(2,2), interpolation="bilinear")(decode_16)
    encode_32 = convolution_block(encode_32, num_filters=48, kernel_size=1)
    decode_32 = layers.Concatenate(axis=-1)([decode_32, encode_32])    
    for i in range(3):
        decode_32 = residual_block(num_filters=256, SE=False)(decode_32)
    
    decode_64 = layers.UpSampling2D(size=(2,2), interpolation="bilinear")(decode_32)    
    decode_64 = residual_block(num_filters=128, SE=False)(decode_64)
    decode_64 = residual_block(num_filters=64, SE=False)(decode_64)
    
    x = layers.UpSampling2D(size=(8,8), interpolation="bilinear")(decode_64)
    
    model_output = layers.Conv2D(num_classes, kernel_size=(1, 1), padding="same", dtype="float32")(x)
    return keras.Model(inputs=model_input, outputs=model_output)