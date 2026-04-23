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


def PrepareTrainValPath(images_dir="Images",
                        masks_dir="Labels_idx",
                        train_ratio=0.8,
                        image_extension="jpg",
                        images_path=None,
                        masks_path=None,
                        seed=1):
    """
    Retrieves file paths for training and validation images.
    - Specify either lists of brightfield image paths and mask image paths via images_path and masks_path,
      or directories containing them via images_dir and masks_dir.
    - The proportion of training data is set with train_ratio; the remainder is used for validation.
    - Paths are randomly sampled, so specify a seed value if reproducibility is required.
    - Specify the file extension of brightfield images with image_extension.
    - Returns four lists: training image paths, training mask paths,
      validation image paths, and validation mask paths.
    """
    if images_path is None:
        ipaths = sorted(glob(f"{images_dir}/*{image_extension}"))
    else:
        ipaths = sorted(images_path)
    
    if masks_path is None:
        mpaths = sorted(glob(f"{masks_dir}/*png"))
    else:
        mpaths = sorted(masks_path)
    
    ipaths = np.array(ipaths)
    mpaths = np.array(mpaths)
    
    random.seed(seed)
    train_ratio = int(len(ipaths)*train_ratio)
    train_idx = random.sample(range(len(ipaths)), train_ratio)
    
    train_ipaths = ipaths[train_idx]
    train_mpaths = mpaths[train_idx]
    val_ipaths = np.delete(ipaths, train_idx)
    val_mpaths = np.delete(mpaths, train_idx)
    
    return train_ipaths, train_mpaths, val_ipaths, val_mpaths
    


def GenerateDataset(
    image_paths:list[str],
    mask_paths:list[str],
    batch_size:int,
    image_size:int,           
    center_zero:bool=True,
    data_augment:bool=True):
    """
    Defines an image generator and creates a TensorFlow dataset.
    - Provide lists of brightfield image paths and mask image paths via image_paths and mask_paths.
    - Data augmentation includes horizontal and vertical flipping, rotation, and translation.
      It is recommended to set data_augment=True when creating the training dataset
      and data_augment=False when creating the validation dataset.
    - Specify the input image size in pixels with image_size.
    - Specify the number of images per mini-batch with batch_size.
    - center_zero: if True, converts pixel values from 0–255 to -1 to 1;
      if False, converts them to 0 to 1.
    """
    def read_image(image_path, image_size=image_size, mask=False, center_zero=True):
        image = tf.io.read_file(image_path)
        if mask:
            image = tf.image.decode_png(image, channels=1)
            image.set_shape([None, None, 1])
            image = tf.image.resize(images=image, size=[image_size, image_size])

        else:
            image = tf.image.decode_png(image)
            image.set_shape([None, None, 3])
            image = tf.image.resize(images=image, size=[image_size, image_size])
            # image = tf.image.random_hue(image=image, max_delta=0.025, seed=None)
            # image = tf.image.random_brightness(image, max_delta=0.05)
            # image = tf.image.random_contrast(image, lower=0.8, upper=1.2)
            # image = tf.image.random_saturation(image,lower=0.8,upper=1.2)
            if center_zero is True:
                image = image / 127.5 - 1
            else:
                image = image/255
        return image
    
    def load_data(image_path, mask_path):
        image = read_image(image_path)
        mask = read_image(mask_path, mask=True)
        return image, mask
    
    def augment_using_layers(image, mask):
        def aug():
            flip = tf.keras.layers.RandomFlip(mode="horizontal_and_vertical")
            rota = tf.keras.layers.RandomRotation(0.2, fill_mode='constant')
            trans = tf.keras.layers.RandomTranslation(height_factor=(-0.2, 0.2),
                                                width_factor=(-0.2, 0.2), 
                                                fill_mode='constant')        
            layers = [flip, trans, rota]
            aug_model = tf.keras.Sequential(layers)
            return aug_model
        aug = aug()
        image_mask = tf.concat([image, mask], -1)  
        image_mask = aug(image_mask)
        image = image_mask[:,:,0:3]
        mask = image_mask[:,:,3]
        return image, tf.cast(mask, 'uint8')
    
    num_imgs = len(image_paths)
    dataset = tf.data.Dataset.from_tensor_slices((image_paths, mask_paths))
    dataset = dataset.map(load_data, num_parallel_calls=tf.data.AUTOTUNE)

    if data_augment is True:
        dataset = dataset.map(lambda x, y: augment_using_layers(x, y), num_parallel_calls=tf.data.AUTOTUNE)
    dataset = dataset.shuffle(buffer_size=num_imgs)
    dataset = dataset.batch(batch_size, drop_remainder=True)
    dataset = dataset.prefetch(buffer_size=tf.data.AUTOTUNE)
    
    return dataset

def GenerateDataset2(
    image_paths:list[str],
    mask_paths:list[str],
    batch_size:int,
    image_size:int,           
    center_zero:bool=True,
    data_augment:bool=True):
    """
    Defines an image generator and creates a TensorFlow dataset.
    - Provide lists of brightfield image paths and mask image paths via image_paths and mask_paths.
    - Data augmentation includes horizontal and vertical flipping, rotation, and translation.
      It is recommended to set data_augment=True when creating the training dataset
      and data_augment=False when creating the validation dataset.
    - Specify the input image size in pixels with image_size.
    - Specify the number of images per mini-batch with batch_size.
    - center_zero: if True, converts pixel values from 0–255 to -1 to 1;
      if False, converts them to 0 to 1.
    """
    def read_image(image_path, image_size=image_size, mask=False, center_zero=True):
        image = tf.io.read_file(image_path)
        if mask:
            image = tf.image.decode_png(image, channels=1)
            image.set_shape([None, None, 1])
            image = tf.image.resize(images=image, size=[image_size, image_size])

        else:
            image = tf.image.decode_png(image)
            image.set_shape([None, None, 3])
            image = tf.image.resize(images=image, size=[image_size, image_size])
            # image = tf.image.random_hue(image=image, max_delta=0.025, seed=None)
            # image = tf.image.random_brightness(image, max_delta=0.05)
            # image = tf.image.random_contrast(image, lower=0.8, upper=1.2)
            # image = tf.image.random_saturation(image,lower=0.8,upper=1.2)
            if center_zero is True:
                image = image / 127.5 - 1
            else:
                image = image/255
        return image
    
    def load_data(image_path, mask_path):
        image = read_image(image_path)
        mask = read_image(mask_path, mask=True)
        return image, mask
    
    def augment_using_layers(image, mask):
        def aug():
            flip = tf.keras.layers.RandomFlip(mode="horizontal_and_vertical")
            rota = tf.keras.layers.RandomRotation(0.2, fill_mode='constant')
            trans = tf.keras.layers.RandomTranslation(height_factor=(-0.2, 0.2),
                                                width_factor=(-0.2, 0.2), 
                                                fill_mode='constant')        
            layers = [flip, trans, rota]
            aug_model = tf.keras.Sequential(layers)
            return aug_model
        aug = aug()
        image_mask = tf.concat([image, mask], -1)  
        image_mask = aug(image_mask)
        image = image_mask[:,:,0:3]
        mask = image_mask[:,:,3]
        return image, tf.cast(mask, 'uint8')
    
    num_imgs = len(image_paths)
    dataset = tf.data.Dataset.from_tensor_slices((image_paths, mask_paths))
    dataset = dataset.map(load_data, num_parallel_calls=tf.data.AUTOTUNE)
    dataset = dataset.cache(filename = "./cache.tf-data")
    if data_augment is True:
        dataset = dataset.map(lambda x, y: augment_using_layers(x, y), num_parallel_calls=tf.data.AUTOTUNE)

    dataset = dataset.shuffle(buffer_size=num_imgs)
    dataset = dataset.batch(batch_size, drop_remainder=True)
    dataset = dataset.prefetch(buffer_size=tf.data.AUTOTUNE)
    
    return dataset

"""
Converts label indices back to an RGB color image.
The colormap is assumed to be an .npy file generated by preparing_tools.MakeColorMap.
To re-import the image as QuPath annotations, an 8-bit color PNG is required.
If palette_mode=True, the image is converted to Pillow palette mode.
"""
def LabelToColorMask(label_img=None,
                     colormap:str=None,
                     label_img_path:str=None,
                     colormap_path:str=None,
                     pallete_mode:bool=True
                     ):
    
    if label_img is None:
        label_img = Image.open(label_img_path)
        label_img = np.asarray(label_img)
    
    rgb = np.zeros(shape=(label_img.shape[0],label_img.shape[1],3), dtype="uint8")
    
    if colormap is None:
        colormap = np.load(colormap_path)
    
    num_classes = len(colormap)
    for i in range(num_classes):
        rgb[label_img==i] = colormap[i]
    
    if pallete_mode is True:
        rgb = Image.fromarray(rgb)
        rgb = rgb.convert(mode="P", matrix=None, dither=0, colors=256, palette=0)
        # rgb.save(f"{save_dir}/{save_name}")

    return rgb

def CheckDataset(dataset, figsize:int=20, colormap_path:str=None):
    """
    Displays samples from the dataset to verify data augmentation.
    - Given a dataset, retrieves 10 batches and displays the first image from each batch.
    - Useful for checking whether data augmentation is functioning as expected.
    """
    plt.figure(figsize=(figsize,figsize))
    i = 1
    for image, mask in dataset.take(10):
        plt.subplot(5, 5, i)
        plt.imshow(keras.preprocessing.image.array_to_img(image[0]))
        plt.grid(False); plt.xticks([]); plt.yticks([]); plt.grid(False)
        i += 1
        plt.subplot(5, 5, i)
        tmp = mask[0]
        if colormap_path is not None:
            tmp = LabelToColorMask(label_img=tmp, colormap_path=colormap_path,pallete_mode=False)
        plt.imshow(tmp)
        plt.grid(False); plt.xticks([]); plt.yticks([]); plt.grid(False)
        i += 1
    plt.show()
    

class UpdatedMeanIoU(tf.keras.metrics.MeanIoU):
    def __init__(self,
                 y_true=None,
                 y_pred=None,
                 num_classes=None,
                 name=None,
                 dtype=None):
        super().__init__(num_classes = num_classes,name=name, dtype=dtype)
    def update_state(self, y_true, y_pred, sample_weight=None):
        y_pred = tf.math.argmax(y_pred, axis=-1)
        return super().update_state(y_true, y_pred, sample_weight)
    
    
def CompileModel(model,
                 num_classes,
                 learning_rate=0.001, 
                 decay=0.0001,
                 IoU=True):
    """
    Compiles the defined model.
    - Specify the model to compile using the model argument.
    - num_classes: number of classification classes
    - Uses Adam optimizer; learning_rate and decay are parameters for Adam.
    - Set IoU=True if IoU is to be used as a metric.

    - Returns the compiled model, so it can be used as:
      model = CompileModel(...)
      (In practice, compile is applied in-place even without reassignment.)
    """
    if IoU==True:
        model.compile(
            optimizer=keras.optimizers.Adam(learning_rate=learning_rate, decay=decay),
            loss=keras.losses.SparseCategoricalCrossentropy(from_logits=True),
            metrics=[UpdatedMeanIoU(num_classes=num_classes, name="MeanIoU"), "accuracy"]
        )
    else:
        model.compile(
            optimizer=keras.optimizers.Adam(learning_rate=learning_rate, decay=decay),
            loss=keras.losses.SparseCategoricalCrossentropy(from_logits=True),
            metrics=["accuracy"]
        )
    return model



def FitModel(
    model,
    train_dataset, 
    val_dataset,
    epochs:int,
    model_save_name:str,
    EarlyStop=True, 
    EarlyStop_patience:int=50, 
    EarlyStop_monitor:str="val_loss", 
    Reduce_lr:bool=True,
    Reduce_lr_monitor:str="val_loss",
    Reduce_lr_factor:float=0.5,
    Reduce_lr_min:float = 1e-12,
    Reduce_lr_patience:int=5,
    Checkpoint:bool=True, 
    Checkpoint_dir:str="training_cp",
    Checkpoint_freq:int=10,
    TensorBoard:bool=True
):
    """
    Trains the model.
    - Specify a compiled model with the model argument.
    - Provide dataset instances created by GenerateDataset() for train_dataset and val_dataset.
    - Set the maximum number of epochs with epochs.
    - model_save_name specifies the filename for saving the trained model.
    - To enable early stopping, set EarlyStop=True.
      Specify the monitored metric with EarlyStop_monitor and the patience (number of epochs with no improvement) with EarlyStop_patience.
    - To adjust the learning rate during training, set Reduce_lr=True.
      Specify the monitored metric with Reduce_lr_monitor and the patience with Reduce_lr_patience.
      Set the reduction factor with Reduce_lr_factor and the minimum learning rate with Reduce_lr_min.
      For example, if Reduce_lr_factor=0.5, the learning rate is halved each time.
    - To enable checkpointing, set Checkpoint=True.
      Specify the save directory and the interval (in epochs) for saving.

    - Returns the training history, so it can be used as:
      history = FitModel(...)
      (Model weights are properly saved during training.)
    """
    callbacks=[]
    if EarlyStop is True:
        earlystop = tf.keras.callbacks.EarlyStopping(monitor=EarlyStop_monitor,patience=EarlyStop_patience)
        callbacks.append(earlystop)
        
    if Reduce_lr is True:
        reduce_lr = keras.callbacks.ReduceLROnPlateau(monitor=Reduce_lr_monitor,factor=Reduce_lr_factor,patience=Reduce_lr_patience,
                                                      verbose=1, mode="auto", min_lr=Reduce_lr_min)
        callbacks.append(reduce_lr)
        
    if Checkpoint is True:
        checkpoint_path = f"{Checkpoint_dir}/cp-{epochs:03d}.ckpt"
        len(list(train_dataset))
        cp_callback = tf.keras.callbacks.ModelCheckpoint(filepath=checkpoint_path, 
                                                         verbose=1, 
                                                         save_weights_only=True,
                                                         save_freq=Checkpoint_freq* len(list(train_dataset)))
        callbacks.append(cp_callback)
    if TensorBoard is True:
        if not os.path.exists("tflog"): os.makedirs("tflog")
        tensorboard_callback = keras.callbacks.TensorBoard(log_dir="tflog/", histogram_freq=1)
        callbacks.append(tensorboard_callback)
    
    if not callbacks:
        history = model.fit(train_dataset,
                        validation_data=val_dataset,
                        epochs=epochs)
    else:
        history = model.fit(train_dataset,
                        validation_data=val_dataset,
                        callbacks=callbacks,
                        epochs=epochs)
    
    model.save(model_save_name)
    model.save(f"{model_save_name}.h5")
    
    return history
    


""" Plot training history """
def PlotHistory(
    history,
    figsize:int=5):
    
    keys = list(history.history.keys())
    keys = [key for key in keys if "lr" not in key]
    ncols = int(len([key for key in keys if "val" not in key]))
    
    fig, axes = plt.subplots(ncols=ncols, figsize=(figsize*ncols,figsize), tight_layout=True, facecolor="whitesmoke")
    
    for i in range(ncols):
        axes[i].plot(history.history[keys[i]]); axes[i].plot(history.history[keys[i+ncols]]); axes[i].set_title(keys[i])
        axes[i].set_ylabel(keys[i]); axes[i].set_xlabel('epoch'); axes[i].legend(['train', 'validation'], loc='upper left')
    
    plt.show()


""" Run prediction on training and validation images """
def CheckPrediction(
    model,
    train_image_paths:list[str],
    val_image_paths:list[str], 
    center_zero:bool=True,
    shuffle:bool=True,
    figsize=5,
    check_num:int=10,
    colormap_path:str=None):
    
    if shuffle:
        random.shuffle(train_image_paths)
        random.shuffle(val_image_paths)
                       
    def read_pred_image(image_path, center_zero=center_zero):
        image = tf.io.read_file(image_path)
        image = tf.image.decode_png(image)
        image = image[:,:,:3]
        image.set_shape([None, None, 3])
        image = tf.cast(image, tf.float32)
        p_img = image/127.5 - 1 if center_zero is True else image/255
        return image, p_img
    
    def show_img_mask(
        path_list:list[str],
        colormap_path=colormap_path,
        figsize=5,
        check_num=10):
        
        images = []; predictions = []
        for i in range(check_num):
            img, p_img = read_pred_image(path_list[i])
            img = keras.preprocessing.image.array_to_img(img)
            images.append(img)
            
            p_img = model.predict(np.expand_dims((p_img), axis=0), verbose=0)
            p_img = np.squeeze(p_img)
            p_img = np.argmax(p_img, axis=-1)
            if colormap_path is not None:
                p_img = LabelToColorMask(label_img=p_img, colormap_path=colormap_path,pallete_mode=False)
                predictions.append(p_img)
            else:
                predictions.append(p_img)
        
        plt.figure(figsize = (figsize,figsize), tight_layout=True)
        j = 1
        for i in range(check_num):
            plt.subplot(check_num,10,j)
            plt.imshow(images[i])
            plt.grid(False); plt.xticks([]); plt.yticks([]); plt.grid(False)
            plt.subplot(check_num,10,j+check_num)
            plt.imshow(predictions[i]); plt.grid(False); plt.xticks([]); plt.yticks([]); plt.grid(False)
            j += 1
        plt.show()

                
    print("Training images")
    show_img_mask(
        path_list=train_image_paths,
        colormap_path=colormap_path,
        figsize=figsize,
        check_num=check_num)    
        
    print("Validation images")
    show_img_mask(
        path_list=val_image_paths,
        colormap_path=colormap_path,
        figsize=figsize,
        check_num=check_num)  

        