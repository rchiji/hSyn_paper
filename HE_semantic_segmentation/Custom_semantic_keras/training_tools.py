## 221221 label画像を表示させる際、元のcolormapに変換して表示できるように。
## 240414 GenerateDataset関数 -> 明視野画像がtileSize以下の場合に、paddingされるように変更

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
    """""""""""""""""""""""
    Train, Val画像のパスを取得
    ・明視野画像のパス一覧、マスク画像のパス一覧をimages_path, masks_path引数に指定するか、
    明視野画像が入っているフォルダへのパス、マスク画像が入っているフォルダへのパスをimages_dir、mask_dir引数に指定
    ・訓練画像の割合をtrain_ratio引数で指定。残りがval画像になる。
    ・画像パスからランダムに取り出すので再現性が必要であればseed引数に任意の値を入れる。
    ・明視野画像の拡張子をimage_extension引数に指定
    ・返り値は訓練明視野画像のパス、訓練マスク画像のパス、検証明視野パス、検証マスク画像パスの4つ
    """""""""""""""""""""""
    ## path一覧の指定が無ければimages_dir, masks_dirからパス一覧を取得
    if images_path is None:
        ipaths = sorted(glob(f"{images_dir}/*{image_extension}"))
    else:
        ipaths = sorted(images_path)
    
    if masks_path is None:
        mpaths = sorted(glob(f"{masks_dir}/*png"))
    else:
        mpaths = sorted(masks_path)
    
    # index操作のためにnumpyに変換
    ipaths = np.array(ipaths)
    mpaths = np.array(mpaths)
    
    # ランダムにindexを取り出す
    random.seed(seed)
    train_ratio = int(len(ipaths)*train_ratio)
    train_idx = random.sample(range(len(ipaths)), train_ratio)
    
    # indexを使って、ファイルパスを分割
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
    """""""""""""""""""""""
    画像生成器を定義してtf用datasetを作るところまで
    ・image_paths, mask_pathsに明視野パス一覧、マスク画像パス一覧を加える。
    ・データ拡張には上下左右反転、回転、平行移動がある。
     訓練dataset作成時にはdata_augment=Trueにし、検証dataset作成時にはdata_augment=Falseにすると良い。
    ・image_sizeにinput画像のpx数を指定
    ・batch_sizeにミニバッチ学習時の１バッチにおける画像数を指定
    ・center_zero: trueなら0-255を-1~1の値に、falseなら0~1の値に変換。
    """""""""""""""""""""""
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
            # keras.layersシリーズのデータ拡張を設定
            flip = tf.keras.layers.RandomFlip(mode="horizontal_and_vertical")
            rota = tf.keras.layers.RandomRotation(0.2, fill_mode='constant')
            trans = tf.keras.layers.RandomTranslation(height_factor=(-0.2, 0.2),
                                                width_factor=(-0.2, 0.2), 
                                                fill_mode='constant')        
            # いったんリストにまとめたものをkeras.Sequentialに入れる。
            layers = [flip, trans, rota]
            aug_model = tf.keras.Sequential(layers)
            return aug_model # <-- keras modelを返す
        aug = aug() # <-- keras modelとして取り出す。
        # 3 chの画像と1 chのマスクを最終軸で重ねる。
        image_mask = tf.concat([image, mask], -1)  
        # データ拡張を実施
        image_mask = aug(image_mask)  # <-- model sequentialを通す。
        # 0-2 chを取り出してimage, 3 chを取り出してmaskとする。
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

## datasetを高速化
# https://zenn.dev/tokyoyoshida/articles/5c3270ce0d4c91#%E7%94%BB%E5%83%8F%E3%82%92%E8%AA%AD%E3%81%BF%E8%BE%BC%E3%82%93%E3%81%A7%E8%A1%A8%E7%A4%BA%E3%81%99%E3%82%8B
def GenerateDataset2(
    image_paths:list[str],
    mask_paths:list[str],
    batch_size:int,
    image_size:int,           
    center_zero:bool=True,
    data_augment:bool=True):
    """""""""""""""""""""""
    画像生成器を定義してtf用datasetを作るところまで
    ・image_paths, mask_pathsに明視野パス一覧、マスク画像パス一覧を加える。
    ・データ拡張には上下左右反転、回転、平行移動がある。
     訓練dataset作成時にはdata_augment=Trueにし、検証dataset作成時にはdata_augment=Falseにすると良い。
    ・image_sizeにinput画像のpx数を指定
    ・batch_sizeにミニバッチ学習時の１バッチにおける画像数を指定
    ・center_zero: trueなら0-255を-1~1の値に、falseなら0~1の値に変換。
    """""""""""""""""""""""
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
            # keras.layersシリーズのデータ拡張を設定
            flip = tf.keras.layers.RandomFlip(mode="horizontal_and_vertical")
            rota = tf.keras.layers.RandomRotation(0.2, fill_mode='constant')
            trans = tf.keras.layers.RandomTranslation(height_factor=(-0.2, 0.2),
                                                width_factor=(-0.2, 0.2), 
                                                fill_mode='constant')        
            # いったんリストにまとめたものをkeras.Sequentialに入れる。
            layers = [flip, trans, rota]
            aug_model = tf.keras.Sequential(layers)
            return aug_model # <-- keras modelを返す
        aug = aug() # <-- keras modelとして取り出す。
        # 3 chの画像と1 chのマスクを最終軸で重ねる。
        image_mask = tf.concat([image, mask], -1)  
        # データ拡張を実施
        image_mask = aug(image_mask)  # <-- model sequentialを通す。
        # 0-2 chを取り出してimage, 3 chを取り出してmaskとする。
        image = image_mask[:,:,0:3]
        mask = image_mask[:,:,3]
        return image, tf.cast(mask, 'uint8')
    
    num_imgs = len(image_paths)
    dataset = tf.data.Dataset.from_tensor_slices((image_paths, mask_paths))
    dataset = dataset.map(load_data, num_parallel_calls=tf.data.AUTOTUNE)
    # ここでローカルにcacheしてみる。一応データ拡張前が良いかなと思って。
    dataset = dataset.cache(filename = "./cache.tf-data")
    if data_augment is True:
        dataset = dataset.map(lambda x, y: augment_using_layers(x, y), num_parallel_calls=tf.data.AUTOTUNE)

    dataset = dataset.shuffle(buffer_size=num_imgs)
    dataset = dataset.batch(batch_size, drop_remainder=True)
    dataset = dataset.prefetch(buffer_size=tf.data.AUTOTUNE)
    
    return dataset

"""
label_indexをRGBのカラー画像に戻す。
colormapはpreparing_tools.MakeColorMapで作成したnpyファイルを想定。
QuPathのアノテーションとして再インポートするためには、8-bit colorのPNGが必要
pallet_mode=TrueでPillowのpalleteモードに変換
"""
def LabelToColorMask(label_img=None,
                     colormap:str=None,
                     label_img_path:str=None,
                     colormap_path:str=None,
                     pallete_mode:bool=True
                     ):
    
    # 画像ファイルの指定が無ければ、画像パスから読み込む
    if label_img is None:
        label_img = Image.open(label_img_path)
        label_img = np.asarray(label_img)
    
    # imgと同じサイズのゼロ行列を用意。ch数は3。
    rgb = np.zeros(shape=(label_img.shape[0],label_img.shape[1],3), dtype="uint8")
    
    # colormapの指定が無ければcolormap_pathから読み込む
    if colormap is None:
        colormap = np.load(colormap_path)
    
    num_classes = len(colormap)
    # 各label index値のところに対応するRGB値を入れていく
    for i in range(num_classes):
        rgb[label_img==i] = colormap[i]
    
    # pallete_mode引数がtrueならrgbを8-bit colorに変換する。
    if pallete_mode is True:
        rgb = Image.fromarray(rgb)
        rgb = rgb.convert(mode="P", matrix=None, dither=0, colors=256, palette=0)
        # rgb.save(f"{save_dir}/{save_name}")

    return rgb

def CheckDataset(dataset, figsize:int=20, colormap_path:str=None):
    """""""""""""""""""""""
    datasetから一部を表示してデータ拡張などを確認
    ・作成したdatasetを指定すると10回分のバッチを取り出して、その1枚目の画像を表示する。
    ・データ拡張が十分か確認できる。
    """""""""""""""""""""""
    plt.figure(figsize=(figsize,figsize))
    i = 1
    for image, mask in dataset.take(10): # data_generatorから10回取り出し。
        plt.subplot(5, 5, i)
        plt.imshow(keras.preprocessing.image.array_to_img(image[0])) # batch単位で取り出されるのでミニバッチの1枚目を表示。
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
    """
    IoUをmonitoringに使う用
    """
    def __init__(self,
                 y_true=None,
                 y_pred=None,
                 num_classes=None,
                 name=None,
                 dtype=None):
        super().__init__(num_classes = num_classes,name=name, dtype=dtype)
    def update_state(self, y_true, y_pred, sample_weight=None):
        y_pred = tf.math.argmax(y_pred, axis=-1)
        # 親クラスに既にあるupdate_state関数を呼び出す。
        # デフォルトではy_predはargmaxされてないものを受け取ることになっているので、上でy_predを処理したものを親クラスのupdate_stateに渡すように工夫。
        return super().update_state(y_true, y_pred, sample_weight)
    
    
def CompileModel(model,
                 num_classes,
                 learning_rate=0.001, 
                 decay=0.0001,
                 IoU=True):
    """""""""""""""""""""""
    定義したモデルをコンパイル
    ・model引数に定義したmodelを指定
    ・分類クラス数: num_classes
    ・optimizerはAdamを使用。learning_rate, decayはAdam用の引数
    ・IoUをmetricsに使用するなら、Trueにする。

    ・返り値はmodelなのでmodel = CompileModel(引数)のように使う。
     (実際はmodel=に返り値を返さなくても引数に指定したモデルのcompileは保存される。)
    """""""""""""""""""""""
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
    """""""""""""""""""""""
    モデルの訓練
    ・model引数にコンパイル済みのmodelを指定
    ・train_dataset, val_datasetにGenerateDataset()で作成したdatasetインスタンスを指定
    ・epochsに最大エポック数を指定
    ・model_save_nameには訓練済みモデルの保存名を記入。
    ・早期終了を仕込む場合、EarlyStop=Trueにする。
     どのmetricsをモニタリング対象とするかをEarlyStop_monitor引数で指定し、monitoring metricsが何epoch変動無ければ早期終了とるつかをEarlyStop_patience引数で指定する。
    ・学習率を変動させる場合、Reduce_lr=Trueにする。
     どのmetricsをモニタリング対象とするかをReduce_lr_monitor引数で指定し、monitoring metricsが何epoch変動が無ければ学習率を変動させるかをReduce_lr_patience引数で指定する。
     学習率の変動係数をReduce_lr_factor, 最小学習率をReduce_lr_minで指定。Reduce_lr_factorが0.5なら学習率が1/2倍ずつ減少する。
    ・学習の途中保存が必要であればCheckpoint=Trueとする。
     保存先フォルダと何epochに一度保存するか指定する。

    ・返り値は学習過程なので、history=FitModel(引数)のように使用する。
     (ちゃんとmodelに重みは保存される。)
    """""""""""""""""""""""
    ## callbackの設定
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
    
    ## callbacksが無かったら。
    if not callbacks:
        history = model.fit(train_dataset,
                        validation_data=val_dataset,
                        epochs=epochs)
    ## callbacksがあったら
    else:
        history = model.fit(train_dataset,
                        validation_data=val_dataset,
                        callbacks=callbacks,
                        epochs=epochs)
    
    ## モデルの保存
    model.save(model_save_name)
    model.save(f"{model_save_name}.h5")
    
    return history
    


""" 学習historyのplot"""
def PlotHistory(
    history,
    figsize:int=5):
    
    keys = list(history.history.keys())
    keys = [key for key in keys if "lr" not in key]
    ncols = int(len([key for key in keys if "val" not in key])) # valとついていないkeyの数を調べる。
    
    fig, axes = plt.subplots(ncols=ncols, figsize=(figsize*ncols,figsize), tight_layout=True, facecolor="whitesmoke")
    
    for i in range(ncols):
        axes[i].plot(history.history[keys[i]]); axes[i].plot(history.history[keys[i+ncols]]); axes[i].set_title(keys[i])
        axes[i].set_ylabel(keys[i]); axes[i].set_xlabel('epoch'); axes[i].legend(['train', 'validation'], loc='upper left')
    
    plt.show()


""" Train, Val画像で予測"""
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
        # 元画像と、予測用にpreprocessした画像を返す。
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
            
            # preprocessした画像を使ってマスクを予測
            p_img = model.predict(np.expand_dims((p_img), axis=0), verbose=0)
            p_img = np.squeeze(p_img)
            p_img = np.argmax(p_img, axis=-1)
            # label indexをRGBに変換。colormap_pathの指定が無ければlabel_indexのまま返す。
            if colormap_path is not None:
                p_img = LabelToColorMask(label_img=p_img, colormap_path=colormap_path,pallete_mode=False)
                predictions.append(p_img)
            else:
                predictions.append(p_img)
        
        # 元画像と予測マスクを並べて表示
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

        