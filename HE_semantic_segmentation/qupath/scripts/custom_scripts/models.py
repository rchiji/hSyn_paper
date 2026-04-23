import os
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'
import tensorflow as tf

def convolution_unit(
    filters:int,   
    kernel_size:int=3,
    strides:int=1,
    name=None,
    dilation_rate=1,
    use_bias=False,
    activation="relu"
):
    if name is not None:
        name_conv,name_bn,name_act = f"{name}_Conv",f"{name}_Bn",f"{name}_Activation"
    else:
        name_conv,name_bn,name_act = None,None,None
        
    def function(inputs):        
        x = tf.keras.layers.Conv2D(
            filters=filters,
            kernel_size=kernel_size,
            strides=strides,
            padding="same",
            kernel_initializer="he_normal",
            use_bias=use_bias,
            dilation_rate=dilation_rate,
            name=name_conv
        )(inputs)
        x = tf.keras.layers.BatchNormalization(name=name_bn)(x)
        if activation is not None:
            x = tf.keras.layers.Activation(activation,name=name_act)(x)
        
        return x
        
    return function


def ResNet50_Unet(
    shape,
    num_classes,
    mixed_float16=True,
    decoder_features = [256,128,64,32,16]
):
    if mixed_float16:
        tf.keras.mixed_precision.set_global_policy('mixed_float16')
        
    model_input = tf.keras.Input(shape=shape)
    
    encoder = tf.keras.applications.ResNet50(
        input_tensor=model_input,
        weights="imagenet",
        include_top=False)
    
    encoder_names = ["conv5_block3_out","conv4_block6_out","conv3_block4_out","conv2_block3_out","conv1_relu"]
    encoder_outputs = [ encoder.get_layer(encoder_name).output for encoder_name in encoder_names]    

    for i in range(len(encoder_outputs)):
        if i == 0:
            x = tf.keras.layers.UpSampling2D(interpolation="bilinear", name=f"Decoder{i}_Upsample")(encoder_outputs[i])
        else:
            x = tf.keras.layers.UpSampling2D(interpolation="bilinear", name=f"Decoder{i}_Upsample")(x)

        if i < len(encoder_outputs)-1:
            x = tf.keras.layers.Concatenate(axis=-1, name=f"Decoder{i}_Concatenate")([x, encoder_outputs[i+1]])
            
        x = convolution_unit(decoder_features[i], name=f"Decoder{i}_ConvUnit1")(x)
        x = convolution_unit(decoder_features[i], name=f"Decoder{i}_ConvUnit2")(x)
    
    x = tf.keras.layers.Conv2D(num_classes, kernel_size=(3,3), padding="same", name="Decoder_final_Conv")(x)
    output = tf.keras.layers.Softmax(axis=-1,dtype=tf.float32, name="Decoder_final_output")(x)
    
    model = tf.keras.Model(inputs=model_input, outputs=output)

    return model


def EffiNetV2S_Unet(
    shape,
    num_classes,
    mixed_float16=True,
    decoder_features = [256,128,64,32,16]
):
    if mixed_float16:
        tf.keras.mixed_precision.set_global_policy('mixed_float16')
        
    model_input = tf.keras.Input(shape=shape)
    
    encoder = tf.keras.applications.EfficientNetV2S(
        weights="imagenet", include_top=False, input_tensor=model_input,
        include_preprocessing = False,
    )
    
    encoder_names = ["top_activation","block6a_expand_activation","block4a_expand_activation","block2d_expand_activation","block1b_project_activation"]
    
    encoder_outputs = [ encoder.get_layer(encoder_name).output for encoder_name in encoder_names]   

    for i in range(len(encoder_outputs)):
        if i == 0:
            x = tf.keras.layers.UpSampling2D(interpolation="bilinear",name=f"Decoder{i}_Upsample")(encoder_outputs[i])
        else:
            x = tf.keras.layers.UpSampling2D(interpolation="bilinear",name=f"Decoder{i}_Upsample")(x)

        if i < len(encoder_outputs)-1:
            x = tf.keras.layers.Concatenate(axis=-1, name=f"Decoder{i}_Concatenate")([x, encoder_outputs[i+1]])
            
        x = convolution_unit(decoder_features[i], name=f"Decoder{i}_ConvUnit1")(x)
        x = convolution_unit(decoder_features[i], name=f"Decoder{i}_ConvUnit2")(x)
    
    x = tf.keras.layers.Conv2D(num_classes, kernel_size=(3,3), padding="same", name="Decoder_final_Conv")(x)
    output = tf.keras.layers.Softmax(axis=-1,dtype=tf.float32, name="Decoder_final_output")(x)
    
    model = tf.keras.Model(inputs=model_input, outputs=output)

    return model


def residual_unit(
    filters:int=256,
    strides:int=1,
    drop_rate:float=0,
    name=None
):

    if name is not None:
        names = [ f"{name}_{suffix}" for suffix in ["ConvUnit1","ConvUnit2","ConvUnit3","DropOut","skip_AvePool","skip_ConvUnit","skip_ConvUnit2","Add","Activation"]]
    else:
        names = [None]*9
        
    def function(inputs):

        x = convolution_unit(filters//4,kernel_size=1,name=names[0])(inputs) 
        x = convolution_unit(filters//4, strides=strides,name=names[1])(x) 
        x = convolution_unit(filters,kernel_size=1,activation=None,name=names[2])(x)
        if drop_rate > 0:
            x = tf.keras.layers.Dropout(rate=drop_rate, noise_shape=(None,1,1,1),name=names[3])(x)

        shortcut = inputs
        
        if strides == 2:
            shortcut = tf.keras.layers.AveragePooling2D(
                pool_size=(2,2),strides=strides,padding="same",name=names[4])(shortcut)
            shortcut = convolution_unit(filters,kernel_size=1,activation=None,name=names[5])(shortcut)
            
        if shortcut.shape[-1] != x.shape[-1]:
            shortcut = convolution_unit(filters,kernel_size=1,activation=None,name=names[6])(shortcut)

        x = tf.keras.layers.Add(name=names[7])([x, shortcut])
        output = tf.keras.layers.Activation("relu",name=names[8])(x)
        
        return output
    
    return function


def residual_block(
    filters:int=256,
    num_loop:int=3,
    strides:int=1,
    drop_rate:float=0,
    name=None
):

    def function(inputs):
        x = inputs
        for i in range(num_loop):
            _name = f"{name}{i+1}" if name is not None else None
                
            if i == 0 and strides == 2:
                x = residual_unit(filters,strides,drop_rate,name=_name)(x)
            else:
                x = residual_unit(filters,strides=1,drop_rate=drop_rate,name=_name)(x)
        return x

    return function


def se_unit(
    filters:int,
    se_denominator:int=8,
    name=None
):
    if name is not None:
        names = [ f"{name}_{suffix}" for suffix in ["AvePool","Reshape","1x1Conv1","1x1Conv2","Multiply"]]
    else:
        names = [None]*5        
    
    def function(inputs):
        x = tf.keras.layers.GlobalAveragePooling2D(name=names[0])(inputs)
        se_shape = (1,1,x.shape.as_list()[-1])
        x = tf.keras.layers.Reshape(se_shape,name=names[1])(x)

        x = tf.keras.layers.Conv2D(
            filters=filters//se_denominator,
            kernel_size=(1,1),
            kernel_initializer="he_normal",
            padding="same",
            activation="relu",
            name=names[2]
        )(x)
        x = tf.keras.layers.Conv2D(
            filters=filters,
            kernel_size=(1,1),
            kernel_initializer="he_normal",
            padding="same",
            activation="sigmoid",
            name=names[3]
        )(x)

        return tf.keras.layers.multiply([inputs, x],name=names[4])

    return function      

def residual_se_unit(
    filters:int=256,
    strides:int=1,
    drop_rate:float=0,
    name=None
):

    if name is not None:
        names = [ f"{name}_{suffix}" for suffix in ["ConvUnit1","ConvUnit2","ConvUnit3","SE","DropOut","skip_AvePool","skip_ConvUnit","skip_ConvUnit2","Add","Activation"]]
    else:
        names = [None]*10
        
    def function(inputs):

        x = convolution_unit(filters//4,kernel_size=1,name=names[0])(inputs) 
        x = convolution_unit(filters//4, strides=strides,name=names[1])(x) 
        x = convolution_unit(filters,kernel_size=1,activation=None,name=names[2])(x)

        x = se_unit(filters,name=names[3])(x)
        
        if drop_rate > 0:
            x = tf.keras.layers.Dropout(rate=drop_rate, noise_shape=(None,1,1,1),name=names[4])(x)

        shortcut = inputs
        
        if strides == 2:
            shortcut = tf.keras.layers.AveragePooling2D(
                pool_size=(2,2),strides=strides,padding="same",name=names[5])(shortcut)
            shortcut = convolution_unit(filters,kernel_size=1,activation=None,name=names[6])(shortcut)
            
        if shortcut.shape[-1] != x.shape[-1]:
            shortcut = convolution_unit(filters,kernel_size=1,activation=None,name=names[7])(shortcut)

        x = tf.keras.layers.Add(name=names[8])([x, shortcut])
        output = tf.keras.layers.Activation("relu",name=names[9])(x)
        
        return output
    
    return function


def residual_se_block(
    filters:int=256,
    num_loop:int=3,
    strides:int=1,
    drop_rate:float=0,
    se:bool=True,
    name=None
):

    def function(inputs):
        x = inputs
        for i in range(num_loop):
            _name = f"{name}{i+1}" if name is not None else None
                
            if i == 0 and strides == 2:
                x = residual_se_unit(filters,strides,drop_rate,name=_name)(x)
            else:
                x = residual_se_unit(filters,strides=1,drop_rate=drop_rate,name=_name)(x)
        return x

    return function