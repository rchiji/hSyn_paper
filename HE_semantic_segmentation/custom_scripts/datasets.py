import os
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'
import numpy as np
from glob import glob
from tqdm.auto import tqdm
import tensorflow as tf
import matplotlib.pyplot as plt
import random
from collections import Counter

def prepare_image_paths(
    image_dir:str,
    extension:str=None,
):
    if extension is None:
        extension = estimate_extension(image_dir)
        print(f"{extension} is used")
        
    image_paths = sorted(glob(f"{image_dir}/*{extension}"))
    if len(image_paths) == 0:
        extension = estimate_extension(image_dir)
        print(f"extension is changed to {extension}")        
        image_paths = sorted(glob(f"{image_dir}/*{extension}"))
        
    return image_paths
    
def estimate_extension(
    image_dir
):
    # 拡張子のリストを取得
    extensions = [os.path.splitext(file)[1] for file in os.listdir(image_dir)]
    # 最頻の拡張子を取得
    most_common_extension = Counter(extensions).most_common(1)[0][0]
    return most_common_extension


def pad_image(image, tile_size):
    # 縦横の不足ピクセル数
    diff_h = tile_size - tf.shape(image)[0]
    diff_w = tile_size - tf.shape(image)[1]
    # 埋めるピクセル数
    padding = [[0, diff_h], [0, diff_w], [0, 0]]
    # 0埋めの実行
    new_image = tf.pad(tensor=image, paddings=padding, mode="CONSTANT", constant_values=0)
    return new_image  


def read_image(file_path, tile_size, mask=False):
    image = tf.io.read_file(file_path)
    image = tf.image.decode_png(image, channels=1) if mask is True else tf.image.decode_png(image, channels=3)
    # 画像サイズがtile_sizeに満たない場合はpadding
    if (tf.shape(image)[0] < tile_size or tf.shape(image)[1] < tile_size):
        image = pad_image(image,tile_size)
    return image


def load_data(image_path, label_path, tile_size):
    image = read_image(image_path,tile_size)
    label = read_image(label_path,tile_size,mask=True)
    return image, label


def shape_augment(image, label):
    # 明視野画像とラベル画像を重ねておく
    combine = tf.concat(values=[image,label],axis=-1)
    # データ拡張内容を定義
    augment_layer = tf.keras.Sequential([
        tf.keras.layers.RandomFlip(mode="horizontal_and_vertical"),
        tf.keras.layers.RandomRotation(factor=0.2, fill_mode="constant", interpolation="nearest"),
        tf.keras.layers.RandomTranslation(height_factor=(-0.2,0.2),width_factor=(-0.2,0.2), fill_mode="constant", interpolation="nearest")
    ])
    # データ拡張を実行
    combine = augment_layer(combine)
    # 明視野画像、ラベル画像を返す
    image = tf.cast(combine[...,0:3], dtype=tf.uint8)
    label = tf.cast(combine[...,3], dtype=tf.uint8)
    label = tf.expand_dims(label, axis=-1)
    return image, label


def rescale_image(image, centering):
    image = tf.cast(image, dtype=tf.float32)
    image = image/127.5-1 if centering is True else image/255.0
    return image


def squeeze_label(label):
    label = tf.squeeze(label)
    label = tf.ensure_shape(label, [None,None])
    return label    

def onehot_label(label, num_classes):
    label = tf.one_hot(indices=label, depth=num_classes)
    return label

def cast_dtype(image, label, label_dtype):
    image = tf.cast(image, dtype=tf.float32)
    label = tf.cast(label, dtype=label_dtype)
    return image, label


def GenerateDataset(
    image_paths:list[str],
    label_paths:list[str],
    batch_size:int,
    tile_size:int,           
    centering:bool=True,
    data_augment:bool=True,
    training:bool=True,
    label_dtype="float32",
    one_hot:bool=False,
    num_classes:int=None
    ):
    
    # 画像のDataset作成
    dataset = tf.data.Dataset.from_tensor_slices(tensors=(image_paths, label_paths)) \
        .map(map_func= lambda image_path, label_path: load_data(image_path,label_path,tile_size), 
             num_parallel_calls=tf.data.AUTOTUNE)
    
    # データ拡張
    if data_augment and training:
        dataset = dataset.map(map_func=shape_augment, num_parallel_calls=tf.data.AUTOTUNE)

    # 1次元の軸を削除
    dataset = dataset.map(map_func= lambda Image, Label: (Image, squeeze_label(Label)),
                          num_parallel_calls=tf.data.AUTOTUNE)  
    # One-hot encoding
    if one_hot:
        dataset = dataset.map(map_func= lambda Image, Label: (Image, onehot_label(Label, num_classes)),
        num_parallel_calls=tf.data.AUTOTUNE)

    # データ型、値の範囲の変更
    dataset = dataset.map(map_func= lambda Image, Label: cast_dtype(Image,Label,label_dtype),
                          num_parallel_calls=tf.data.AUTOTUNE) \
                     .map(map_func= lambda Image, Label: (rescale_image(Image, centering), Label),
                          num_parallel_calls=tf.data.AUTOTUNE)
	
    if training:
	    dataset = dataset.cache() \
            .shuffle(buffer_size=len(dataset)) \
            .batch(batch_size=batch_size, drop_remainder=True) \
            .prefetch(buffer_size=tf.data.AUTOTUNE)
    else:
        dataset = dataset.batch(batch_size=batch_size, drop_remainder=False) \
                         .prefetch(buffer_size=tf.data.AUTOTUNE)

    return dataset


def CheckDataset(dataset, num:int=10):
    import math
    nrows = math.ceil( num/ 5 ) * 2 # plotの行数設定
    ncols = min(num, 5)
    fig, axs = plt.subplots(ncols=ncols, nrows=nrows)
    i = 0
    mini_dataset = dataset.take(num)
    
    ## One hotかどうか判定
    _,la = next(iter(dataset))
    is_onehot = tf.reduce_all( # 全ての要素が0か1かを判定
        tf.logical_or(tf.equal(la, 0), tf.equal(la, 1))
    )
    
    for images, labels in mini_dataset: 
        row = i // 5*2
        col = i % 5
        
        # 明視野画像
        axs[row,col].imshow(tf.keras.preprocessing.image.array_to_img(images[0]))
        axs[row,col].grid(False); axs[row,col].set_xticks([]); axs[row,col].set_yticks([])
        
        # ラベル画像
        label = labels[0]
        # one-hotの場合は、argmax
        if is_onehot:
            label = tf.argmax(label,axis=-1)
        axs[row+1,col].imshow(label)
        axs[row+1,col].grid(False); axs[row+1,col].set_xticks([]); axs[row+1,col].set_yticks([])
        i += 1
        
    plt.tight_layout()
    plt.show()


def PrepareTrainValPath(image_dir="Images",
                        label_dir="Labels",
                        train_ratio=0.8,
                        image_extension="jpg",
                        label_extension="png",
                        image_paths=None,
                        label_paths=None,
                        seed=None):
    """Train, Val画像のパスを取得する関数
    
    Args:
        image_dir (str): 明視野画像が入っているフォルダのパス
        label_dir (str): ラベル画像が入っているフォルダのパス
        train_ratio (float): 学習画像の割合
        image_extension (str): 明視野画像の拡張子
        label_extension (str): ラベル画像の拡張子
        image_paths (list): 明視野画像のパス一覧
        label_paths (list): ラベル画像のパス一覧
        seed (int): ランダムシード

    Returns:
        tuple: 訓練明視野画像のパス、訓練マスク画像のパス、検証明視野パス、検証マスク画像パスの4つ
    """
    ## path一覧の指定が無ければimage_dir, label_dirからパス一覧を取得
    if image_paths is None:
        ipaths = prepare_image_paths(image_dir, image_extension)
    else:
        ipaths = sorted(image_paths)
    
    if label_paths is None:
        lpaths = prepare_image_paths(label_dir, label_extension)
    else:
        lpaths = sorted(label_paths)
    
    # index操作のためにnumpyに変換
    ipaths = np.array(ipaths)
    lpaths = np.array(lpaths)
    
    # ランダムにindexを取り出す
    if seed is not None:
        random.seed(seed)
    train_ratio = int(len(ipaths)*train_ratio)
    train_idx = random.sample(range(len(ipaths)), train_ratio)
    
    # indexを使って、ファイルパスを分割
    train_ipaths = ipaths[train_idx]
    train_lpaths = lpaths[train_idx]
    val_ipaths = np.delete(ipaths, train_idx)
    val_lpaths = np.delete(lpaths, train_idx)
    
    return train_ipaths, train_lpaths, val_ipaths, val_lpaths

def GenerateDatasets(
    image_dir:str,
    label_dir:str,
    batch_size:int,
    tile_size:int,
    centering=True,
    data_augment=True,
    label_dtype="float32",
    one_hot:bool=False,
    num_classes:int=None,
    train_ratio=0.8,
    image_extension="jpg",
    label_extension="png",
    seed=None
):
    train_image_paths, train_label_paths, val_image_paths, val_label_paths = PrepareTrainValPath(
        image_dir, label_dir, train_ratio, image_extension, label_extension, seed=seed)

    arguments = {"batch_size":batch_size,
                 "tile_size":tile_size, 
                 "centering":centering,
                 "data_augment":data_augment, 
                 "label_dtype":label_dtype,
                 "one_hot":one_hot,
                 "num_classes":num_classes}
    
    train_dataset  = GenerateDataset(image_paths=train_image_paths,
                                     label_paths=train_label_paths,
                                     training=True,
                                     **arguments)

    val_dataset  = GenerateDataset(image_paths=val_image_paths, 
                                   label_paths=val_label_paths,
                                   training=False, 
                                   **arguments)
    
    return train_dataset, val_dataset


def read_pred_image(path):
    """
    ファイルパスから明視野画像を読み込んでtf.float32型で返す
    """
    image = tf.io.read_file(path)
    image = tf.image.decode_png(image, channels=3)
    image = tf.cast(image, tf.float32)
    return image
    
def GeneratePredictDataset(
    image_paths=None,
    image_dir:str=None,
    extension:str="jpg",
    batch_size:int=50,
    tile_size:int=512,           
    centering:bool=True,
    ):
    """
    予測時に使用する明視野画像のdatasetを返す機能
    """

    # 画像パスリストの用意が無ければフォルダからパス作成
    if image_paths is None:
        image_paths = prepare_image_paths(image_dir,extension)
    
    image_paths = sorted(image_paths)
    
    # 明視野画像のDataset作成
    dataset = tf.data.Dataset.from_tensor_slices(tensors=image_paths) \
        .map(lambda path : (read_pred_image(path), path), num_parallel_calls=tf.data.AUTOTUNE) \
        .map(lambda image, path : (pad_image(image, tile_size), path), 
             num_parallel_calls=tf.data.AUTOTUNE) \
        .map(lambda image, path : (rescale_image(image, centering), path),
             num_parallel_calls=tf.data.AUTOTUNE) \
        .batch(batch_size=batch_size, drop_remainder=False) \
        .prefetch(buffer_size=tf.data.AUTOTUNE)

    return dataset


def GeneratePredictDataset_labelholder(
    labelholders,
    batch_size:int=50,
    tile_size:int=512,           
    centering:bool=True,
    ):

    # labelholdersから画像ファイルパスのリスト作成
    image_paths = [ lh.path for lh in labelholders]
    
    # 明視野画像のDataset作成
    dataset = tf.data.Dataset.from_tensor_slices(tensors=image_paths) \
        .map(read_pred_image, num_parallel_calls=tf.data.AUTOTUNE) \
        .map(lambda image : pad_image(image, tile_size), 
             num_parallel_calls=tf.data.AUTOTUNE) \
        .map(lambda image : rescale_image(image, centering),
             num_parallel_calls=tf.data.AUTOTUNE) \
        .batch(batch_size=batch_size, drop_remainder=False) \
        .prefetch(buffer_size=tf.data.AUTOTUNE)

    return dataset


def _add_sample_weights(image, label, class_weight):
    sample_weights = tf.gather(class_weight, indices=tf.cast(label, tf.int32))
    return image, label, sample_weights

def add_sample_weights(dataset, class_weight):
    dataset = dataset.map(lambda Image, Label: _add_sample_weights(Image,Label,class_weight))
    return dataset


def read_region(slide,location,size):
    image = slide.read_region(location=location,
                              level=0,
                              size=size)
    image = tf.keras.utils.img_to_array(image.convert("RGB"))
    image = tf.convert_to_tensor(np.array(image))
    image = tf.cast(image, tf.float32)
    
    return image