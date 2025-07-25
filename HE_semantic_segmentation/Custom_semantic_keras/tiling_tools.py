## 221125 updated

from tqdm.auto import tqdm
import re
from glob import glob
import numpy as np
import cv2
import scipy.ndimage
from numpy import int8
import os
from PIL import Image
Image.MAX_IMAGE_PIXELS = 1000000000

""" == GetSlideNames ==
スライド名一覧を取り出す
デフォルトでは予測マスク画像に対して実行するように、prefixにLabels_が付いた画像を探すようにしている。
Labels_slidename [d=,x=,y=,w=,h=].pngというファイル名を想定。
半角スペースで区切ってslidenameのところを取り出す。
windowsパスに入る\\で画像パスの文字列をsplitするようにしているため, \\以外のディレクトリ区切りの場合ヒットが難しい。
"""
def GetSlideNames(
    path_list:list[str],
    file_dir:str=None,
    prefix:str="Labels_",
    extension:str="png"
):
    
    # 画像のパス一覧の指定がNoneの場合、dir=引数の指定先から画像のパス一覧を取得する。
    if path_list is None:
        path_list = sorted(glob(f"{file_dir}/*{extension}"))
        print(f"Got mask path from {file_dir}")
    else:
        path_list = sorted(path_list)
    
    print("Getting slide names")
    slidenames = []
    for path in path_list:
        path = os.path.basename(path) # ファイル名を取得
        if prefix is None:
            slidenames.append(re.split("[ ]", path)[0]) # 半角スペースで区切ってスライド名と座標情報を分ける。
        else:
            if path.startswith(prefix):
                path = path.replace(prefix, "")
                slidenames.append(re.split("[ ]", path)[0])
        slidenames = sorted(list(set(slidenames)))
    print(f"{len(slidenames)} slides were found")
    return slidenames


""" == MergeTiles ==
1スライドずつタイル画像を張り付けて1枚絵にする 
・path_listに画像パスの一覧を指定するか、画像が入っているディレクトリを指定しても良い。
・タイル切り出し前の画像での座標情報をファイル名から読み取るので[d=,x=,y=,w=,h=]の形式でファイル名に入っている必要がある。
・downsample係数は指定しなくてもd=の箇所から読み取るが、指定してもよい。
(QuPathでは等倍で書き出すとファイル名にd=が残らない。d=が無い場合は等倍と読み取る)
・TrimOverlap=Trueにするとoverlapを考慮して画像の外周を削って中心のみにする
　TrimOverlap=Falseならoverlap_px引数を指定していても0に上書きする。
 (overlap領域の確率マップを平均化しているならtrimしてもしなくても変わらないはず。) 
・保存先フォルダ引数save_dirは存在しなければフォルダが作られる。
・makeColorMaskをTrueにすると、label index画像を元のカラー画像に変換できる。
　その場合はpreparing_tools.MakeColorMapで作成したColorMapのnpyファイルへのパスを入れる。
  label index画像とは異なるフォルダに別で保存するようにしている。
  8bit colorで保存したければpallete_modeをTrueにする。(なぜか背景pxの変換が失敗することがある。)
・npy形式で保存したければsaved_format_npy=Trueにする。Falseなら.pngで保存される。
"""
def MergeTiles(pred_mask_dir:str="Pred_masks",
               imageSize:int=512,
               TrimOverlap:bool=False,
               overlap_px:int=256,
               downsample:int=None,
               save_dir:str="Pred_masks_merge",
               prefix:str="Labels",
               makeColorMask:bool=True,
               save_ColorMask_dir:str="Pred_ColorMasks",
               colormap_path:str="index.npy",
               pallete_mode:bool=False,
               path_list:list[str]=None,
               saved_format_npy:bool=False):
    
    # 予測マスクのパス一覧の指定が無ければpred_mask_dirからnpyファイルパス一覧を取得
    if path_list is None:
        if saved_format_npy is True:
            path_list = sorted(glob(f"{pred_mask_dir}/{prefix}*npy"))
        else:
            path_list = sorted(glob(f"{pred_mask_dir}/{prefix}*png"))
        print(f"Got pred_mask path from {pred_mask_dir}")
    
    ## 保存先フォルダの作成
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    if save_ColorMask_dir is not None:
        if not os.path.exists(save_ColorMask_dir):
            os.makedirs(save_ColorMask_dir)
            
    ## 画像中心のCropナシを指定していればoverlap_px引数は0にする。
    if TrimOverlap is False:
        overlap_px = 0
    
    # 上下両サイドの削るpixel数を定義
    overlap = int(overlap_px/2)
    
    slidenames = sorted(GetSlideNames(path_list=path_list, file_dir=None, extension="png"))
    # スライド名ごとにマージしたマスク画像を作成
    for slide in tqdm(slidenames):
        print(f"<<<<<< Merging {slide} >>>>>>")
        
        # そのスライドのファイルパス一覧を取得
        paths = [s for s in path_list if slide in s]
        
        # ファイル名からダウンサンプリング係数とタイルサイズを保存。
        if downsample is None:
            path = os.path.basename(paths[0])
            # ダウンサンプリング係数が入っていなかったら1
            if not "d=" in path:
                downsample = int(1)
            # ある場合はファイル名のd=のところから読み取る。
            else:
                tile_info = [s for s in re.split("[ ]", path) if s.startswith('[')][-1] # 半角スペースでファイルパスをsplitして[で始まる最後の要素を取り出す。
                tile_info = re.split("[,.\\[\\]]", tile_info) # , . [ ]でsplitしてリスト化
                downsample = int(re.split("[=]", tile_info[1])[1]) # =でsplitして1番目の要素をダウンサンプリング係数として保存 
                
        # w=かh=のところからタイルサイズを読み取る。(元からd=が有るか無いかでずれる)
        if imageSize is None:
            path = os.path.basename(paths[0])
            tile_info = [s for s in re.split("[ ]", path) if s.startswith('[')][-1]
            tile_info = re.split("[,.\\[\\]]", tile_info)
            # ダウンサンプリング係数が掛け算された後のpx sizeがファイル名に入っている。(実際の画像pxサイズではなくダウンサンプルしていない場合のpxサイズ)
            tile_size = int(re.split("[=]", tile_info[4])[1])
        else:
            tile_size = int(imageSize*downsample)

        # ファイル名からx座標, y座標のリストを作る。(ダウンサンプリング前のWSIでの座標がファイル名に初めから入っている。)
        x=[]; y=[]
        for path in paths:    
            x.append(int([p for p in re.split("[,\\[\\]]", os.path.basename(path)) if p.startswith("x")][0][2:])) # x=__で始まる箇所を抜き出して3文字目以降を取得
            y.append(int([p for p in re.split("[,\\[\\]]", os.path.basename(path)) if p.startswith("y")][0][2:]))
        # 最大値を取得
        x_max = max(x); y_max = max(y)

        # タイル結合後のサイズ分の空のnumpyアレイを作成
        merge = np.zeros((y_max + tile_size, x_max + tile_size), dtype = int8)

        # タイル画像を結合していく
        for path in tqdm(paths):
            # タイル画像を読み込み。※ この段階ではダウンサンプリング後のpx width, height
            if saved_format_npy is True:
                img = np.load(path)
            else:
                img = Image.open(path)
                img = np.asarray(img)
                
            # Upsamplingしてダウンサンプリング前の形状に戻す。
            img = scipy.ndimage.zoom(img, downsample, order=0)

            # 上下左右のoverlap分を削って中心をcropする。
            crop_end = int(tile_size- overlap)
            img = img[overlap:crop_end, overlap:crop_end]

            # 画像のファイル名からx, y座標を取り出し。削ったpixel数/2を足して貼り付け座標を調整
            x_start = int([p for p in re.split("[,\\[\\]]", os.path.basename(path)) if p.startswith("x")][0][2:]) + overlap
            y_start = int([p for p in re.split("[,\\[\\]]", os.path.basename(path)) if p.startswith("y")][0][2:]) + overlap
            # xy座標にタイル画像サイズを足して、タイル画像の範囲を取り出す。
            x_end = x_start + tile_size-overlap_px
            y_end = y_start + tile_size-overlap_px
            # タイル画像範囲にタイル画像の輝度値を代入する
            merge[y_start:y_end, x_start:x_end] = img

        # 結合後の画像を再度ダウンサンプリング
        merge = scipy.ndimage.zoom(merge, 1/downsample, order=0)
        
        # QuPathの書き出し時の名前っぽく保存名をつける。
        save_name = f"{slide} [d={downsample},x=0,y=0,w={x_max+(imageSize*downsample)},h={y_max+(imageSize*downsample)}].png"
        # label index画像を保存
        cv2.imwrite(f"{save_dir}/{save_name}", merge)
        
        # label index画像をカラー画像に変換
        if makeColorMask is True:
            print("Converting RGB mask")
            colormask = LabelToColorMask(label_img=merge,
                                         colormap_path=colormap_path,
                                         pallete_mode=pallete_mode)
            if pallete_mode is True:
                colormask.save(f"{save_ColorMask_dir}/{save_name}")
            else:
                # cv2で保存するならBGR2RGBをしてないと色がおかしくなる。grayscaleなら関係ないみたい。
                colormask = cv2.cvtColor(colormask, cv2.COLOR_BGR2RGB)
                cv2.imwrite(f"{save_ColorMask_dir}/{save_name}", colormask)

        print(slide + " is saved!!")

""" == MergeTiles_from_ImageHolder ==
1スライドずつタイル画像を張り付けて1枚絵にする 
predict_tools.pyで定義したImageHolderクラスを使用する版。
座標情報などはImageHolderインスタンスのメンバ変数に記録してあるので、それらの情報を使用する。
"""
def MergeTiles_from_ImageHolder(
    imgholder_list,
    imageSize:int=512,
    TrimOverlap:bool=False,
    overlap_px:int=256,
    downsample:int=None,
    save_dir:str="Pred_masks_merge",
    prefix:str="Labels",
    makeColorMask:bool=True,
    save_ColorMask_dir:str="Pred_ColorMasks",
    colormap_path:str="index.npy",
    pallete_mode:bool=False,
    path_list:list[str]=None,
):
    
    ## 保存先フォルダの作成
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    if save_ColorMask_dir is not None:
        if not os.path.exists(save_ColorMask_dir):
            os.makedirs(save_ColorMask_dir)
    
    slidename = imgholder_list[0].slidename
    
    ## 画像中心のCropナシを指定していればoverlap_px引数は0にする。
    if TrimOverlap is False:
        overlap_px = 0
    # 上下両サイドの削るpixel数を定義
    overlap = int(overlap_px/2)
        
    # ImageHolderクラスのメンバ変数からダウンサンプリング係数とタイルサイズを保存。
    if downsample is None:
        downsample = imgholder_list[0].tile_info["d"]
    tile_size = imgholder_list[0].tile_size

    # ファイル名からx座標, y座標のリストを作る。(ダウンサンプリング前のWSIでの座標がファイル名に初めから入っている。)
    x=[]; y=[]
    for imgholder in imgholder_list:    
        x.append(imgholder.tile_info["x"]); y.append(imgholder.tile_info["y"])
    # 最大値を取得
    x_max = max(x); y_max = max(y)

    # タイル結合後のサイズ分の空のnumpyアレイを作成
    merge = np.zeros((y_max + tile_size, x_max + tile_size), dtype = np.int8)

    # タイル画像を結合していく
    # マスクを結合したらImageHolderインスタンスはもう要らないので、popしてCPUメモリから削除していく
    print("--- Tiling masks ---")
    pbar = tqdm(total=len(imgholder_list)) # while loop用のtqdm設定
    while imgholder_list:
        imgholder = imgholder_list.pop()
        
        img = imgholder.mask_final
        
        # Upsamplingしてダウンサンプリング前の形状に戻す。
        img = scipy.ndimage.zoom(img, downsample, order=0)

        # 上下左右のoverlap分を削って中心をcropする。
        crop_end = int(tile_size - overlap)
        img = img[overlap:crop_end, overlap:crop_end]

        # 画像のファイル名からx, y座標を取り出し。削ったpixel数/2を足して貼り付け座標を調整
        x_start = imgholder.tile_info["x"] + overlap
        y_start = imgholder.tile_info["y"] + overlap
        # xy座標にタイル画像サイズを足して、タイル画像の範囲を取り出す。
        x_end = x_start + tile_size - overlap_px
        y_end = y_start + tile_size - overlap_px
        # タイル画像範囲にタイル画像の輝度値を代入する
        merge[y_start:y_end, x_start:x_end] = img 
        
        pbar.update(1) # while loop用のtqdm設定
    pbar.close() # while loop用のtqdm設定

    # 結合後の画像を再度ダウンサンプリング
    merge = scipy.ndimage.zoom(merge, 1/downsample, order=0)

    # QuPathの書き出し時の名前っぽく保存名をつける。
    save_name = f"{slidename} [d={downsample},x=0,y=0,w={x_max+(imageSize*downsample)},h={y_max+(imageSize*downsample)}].png"
    # label index画像を保存
    cv2.imwrite(f"{save_dir}/{save_name}", merge)

    # label index画像をカラー画像に変換
    if makeColorMask is True:
        print("--- Converting RGB mask ---")
        colormask = LabelToColorMask(label_img=merge,
                                     colormap_path=colormap_path,
                                     pallete_mode=pallete_mode)
        if pallete_mode is True:
            colormask.save(f"{save_ColorMask_dir}/{save_name}")
        else:
            # cv2で保存するならBGR2RGBをしてないと色がおかしくなる。grayscaleなら関係ないみたい。
            colormask = cv2.cvtColor(colormask, cv2.COLOR_BGR2RGB)
            cv2.imwrite(f"{save_ColorMask_dir}/{save_name}", colormask)
        
        
""" == ReSplit ==
タイリングしたマスク画像を再タイル化 (QuPathで扱う際に軽量化させるため)
・ファイル名から元画像での座標情報を読み取るので[d=,x=,y=,w=,h=]がファイル名に入っているようにする。
・downsampleは指定が無ければファイル名から読み取る。
・path_list引数にmerge後の画像パス一覧をpath_list引数に指定してもいいし、merge画像が入っているディレクトリを指定merged_mask_dir引数にしてもよい
・tileSizeは任意。実際はtileSize * downsampleが書き出し画像のpx数になる。
・保存先フォルダ引数save_dirは存在しなければフォルダが作られる。
"""

def ReSplit(merged_mask_dir:str="Pred_masks_merge",
            path_list:list[str]=None,
            tileSize:int=2000, 
            save_dir:str="Pred_masks_split", 
            downsample:int=None):
    
    # 画像パス一覧が無ければmerged_mask_dirから画像一覧パスを取得
    if path_list is None:
        maskPaths = sorted(glob(f"{merged_mask_dir}/*png"))
    else:
        maskPaths = sorted(path_list)
    
    ## 保存先フォルダの作成
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    # スライド名
    for path in tqdm(maskPaths):
        basename = os.path.basename(path)
        slide_name = [s for s in re.split("[ ]", basename)][0]
        print(f"<<< Spliting -- {slide_name} -- >>>")
        
        # ファイル名からダウンサンプリング係数とタイルサイズを保存。
        if downsample is None:
            # ダウンサンプリング係数が入っていなかったら1を
            if not "d=" in basename:
                downsample = int(1)
            # ある場合はファイル名のd=のところから読み取る。
            else:
                tile_info = [s for s in re.split("[ ]", basename) if s.startswith('[')][-1] # 半角スペースでファイルパスをsplitして[で始まる最後の要素を取り出す。
                tile_info = re.split("[,.\\[\\]]", tile_info) # , . [ ]でsplitしてリスト化
                downsample = int(re.split("[=]", tile_info[1])[1]) # =でsplitして1番目の要素をダウンサンプリング係数として保存     

        # 画像読み込み
        mask = np.asarray(Image.open(path))

        # 画像サイズを取得
        height = mask.shape[0]; width = mask.shape[1]

        # 分割時の座標リスト作成
        heightSplit = np.arange(0, height, tileSize)
        heightSplit = np.append(heightSplit , (heightSplit[-1] + height % tileSize)) # 剰余はappendで追加
        widthSplit = np.arange(0, width, tileSize)
        widthSplit = np.append(widthSplit , (widthSplit[-1] + width % tileSize))

        for i, w in enumerate(widthSplit[:-1]):
            for j, h in enumerate(heightSplit[:-1]):
                maskTile = mask[h:heightSplit[j+1], w:widthSplit[i+1]]
                tileHeight = maskTile.shape[0]
                tileWidth = maskTile.shape[1]

                # QuPathの書き出し時の名前っぽく保存
                cv2.imwrite(f"{save_dir}/{slide_name} [d={downsample},x={w},y={h},w={tileWidth*downsample},h={tileHeight*downsample}].png", maskTile)

        print(f"{slide_name} is saved!!")

        
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

## LabelToColorMaskを指定フォルダの画像に一括処理する版
def LabelToColorMask_Batch(
    label_dir:str,
    colormap=None,
    colormap_path:str="index.npy",
    save_dir:str=None,
    pallete_mode:bool=False):
    
    ## 保存先フォルダの作成
    if not os.path.exists(save_dir):
        os.makedirs(save_dir)
    
    # colormapの指定が無ければcolormap_pathから読み込む
    if colormap==None:
        colormap = np.load(colormap_path)
    num_classes = len(colormap)
    print("===== Color map =====")
    print(colormap); print("====================")
    
    imgPath = sorted(glob(f"{label_dir}/*png"))
    print("Converting index to color")
    for path in tqdm(imgPath):
        fname = os.path.basename(path)
        print(fname)
        label_img = Image.open(path)
        label_img = np.asarray(label_img)
    
        # imgと同じサイズのゼロ行列を用意。ch数は3。
        rgb = np.zeros(shape=(label_img.shape[0],label_img.shape[1],3), dtype="uint8")

        # 各label index値のところに対応するRGB値を入れていく
        
        for i in tqdm(range(num_classes)):
            rgb[label_img==i] = colormap[i]
            
        # pallete_mode引数がtrueならrgbを8-bit colorに変換して保存する。Falseなら24-bit colorとして保存。
        if pallete_mode is True:
            rgb = Image.fromarray(rgb)
            rgb = rgb.convert(mode="P", matrix=None, dither=0, colors=256, palette=0)
            rgb.save(f"{save_dir}/{fname}")
        else:
            rgb = cv2.cvtColor(rgb, cv2.COLOR_BGR2RGB)
            cv2.imwrite(f"{save_dir}/{fname}", rgb)

