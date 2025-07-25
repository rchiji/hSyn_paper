import os
os.environ['TF_FORCE_GPU_ALLOW_GROWTH'] = 'true'
import numpy as np
from glob import glob
from tqdm.auto import tqdm
import tensorflow as tf
import matplotlib.pyplot as plt
import math


def PlotHistory(
    history,
    figsize:int=4,
    show_best=True,
    log_base=None,
    epoch_from=None,
    epoch_to=None,
):
    import tensorflow as tf
    import matplotlib.pyplot as plt
    import math
    import copy

    # dict型を取り出しておく
    hist = copy.copy(history)
    if isinstance(hist, tf.keras.callbacks.History):
        hist = hist.history

    # epochの表示域分のデータを抽出
    if epoch_from is not None:
        hist = { key:values[epoch_from:] for key,values in hist.items()}
    if epoch_to is not None:
        hist = { key:values[epoch_to:] for key,values in hist.items()}
                      
    # plotのkeyを取得
    keys = list(hist.keys())
    plot_keys = [key for key in keys if "val" not in key]
    num_plots = int(len(plot_keys))

    ## 適切なfigure layoutを探索
    # 初期値
    nrows, ncols = 1, num_plots
    best_layout = (nrows, ncols)
    # 最も余白が少なくなる、行列数を決定
    if num_plots > 5:
        min_wasted_space = float('inf')
        # nrowとncolの候補を計算
        for nrow in range(2, int(math.sqrt(num_plots)) + 3):
            for max_ncol in [4,5]:
                ncol = min(max_ncol, math.ceil(num_plots / nrow))
                total_subplots = nrow * ncol
                if total_subplots < num_plots:
                    continue
                wasted_space = total_subplots - num_plots
        
                if wasted_space < min_wasted_space:
                    min_wasted_space = wasted_space
                    best_layout = (nrow, ncol)
        nrows, ncols = best_layout
    
    fig, axes = plt.subplots(ncols=ncols, nrows=nrows,squeeze=False,
                             figsize=(figsize*ncols,figsize*nrows),
                             tight_layout=True, facecolor="whitesmoke")

    ## historyをplot
    for i, ax in enumerate(axes.flatten()):
        if i < num_plots:
            name = plot_keys[i]
            # lossを対数化する場合は値を変換
            if "loss" in name and log_base is not None:
                log_values = [math.log(value,log_base) for value in hist[name]]
                train_line, = ax.plot(log_values, label="train")
            # それ以外は元の値をplot
            else:
                train_line, = ax.plot(hist[name], label="train")
            ax.set_title(name)
            ax.set_ylabel(name)
            ax.set_xlabel('epoch')

            # ベストな値と横線の表示
            if show_best:
                # 最大値とインデックスを取得
                if any(term in name for term in ["loss","mse","msge"]):
                    best_value = min(hist[name])
                else:
                    best_value = max(hist[name])
                best_idx = hist[name].index(best_value)
                # 横線を引く
                ax.axhline(y=best_value, color='r', linestyle='--') 
                # 箇所表示
                ax.text(best_idx, best_value, str(round(best_value,3)), ha='center')
            
            val_name = "val_" + name
            if val_name in keys:
                # lossを対数化する場合は値を変換
                if "loss" in val_name and log_base is not None:
                    log_values = [math.log(value,log_base) for value in hist[val_name]]
                    val_line, = ax.plot(log_values,label="val")
                # それ以外は元の値をplot
                else:
                    val_line, = ax.plot(hist[val_name],label="val")
                    
                if show_best:
                    # 最大値とインデックスを取得
                    if any(term in val_name for term in ["loss","mse","msge"]):
                        best_value = min(hist[val_name])
                    else:
                        best_value = max(hist[val_name])
                    best_idx = hist[val_name].index(best_value)
                    # 横線を引く
                    ax.axhline(y=best_value, color='r', linestyle='-.') 
                    # 箇所表示
                    ax.text(best_idx, best_value, str(round(best_value,3)), ha='center') 
                
                ax.legend(handles=[train_line, val_line],
                    loc='upper left')
        else:
            # 余分なサブプロットは非表示に
            ax.axis('off')
    plt.show()


def CheckPrediction(
    model,
    dataset,
    cmap="viridis",
    num_classes=None,
    figsize=3,
    fontsize=20,
    check_num:int=5):
    
    import tensorflow as tf
    import matplotlib.pyplot as plt
    
    if num_classes is not None:
        import matplotlib.colors as colors
        cm = plt.get_cmap(name=cmap, lut=num_classes-1) 
        cmap = colors.ListedColormap([(1.0, 1.0, 1.0)] + [ cm(i)[:3] for i in range(cm.N)])

    # datasetからcheck_num組の画像とラベルを取り出す
    dataset = dataset.unbatch().shuffle(buffer_size=1000)
    dataset = dataset.take(check_num)
    image_list = []; orig_image_list = []; truth_list = []
    is_onehot = None
    for image, label in dataset:
        image_list.append(image)
        
        # one-hot対策
        if is_onehot is None:
            is_onehot = np.all((label.numpy() == 0) | (label.numpy() == 1))     
      
        if is_onehot:
            label = tf.argmax(label,axis=-1)
            
        truth_list.append(label)
        # 値の範囲を元に戻した明視野画像
        orig_image_list.append(tf.keras.preprocessing.image.array_to_img(image))

    # 予測用に画像リストを重ねてバッチ軸を作成
    images = tf.stack(image_list, axis=0)
    # 予測
    preds = model.predict(images)
    # 最も確率が高いラベルへ
    preds = tf.argmax(preds, axis=-1)

    # 元画像/予測画像/正解画像 を並べて表示
    fig, axes = plt.subplots(ncols=check_num, nrows=3,
                             figsize=(figsize*check_num,figsize*3),
                             tight_layout=True)
    
    for i in range(check_num):
        axes[0,i].imshow(orig_image_list[i])
        axes[0,i].grid(False); axes[0,i].set_xticklabels([]); axes[0,i].set_yticklabels([]); plt.grid(False)
        axes[1,i].imshow(truth_list[i], vmin=0, vmax=num_classes, cmap=cmap)
        axes[1,i].grid(False); axes[2,i].set_xticklabels([]); axes[2,i].set_yticklabels([]); plt.grid(False)        
        axes[2,i].imshow(preds[i], vmin=0, vmax=num_classes, cmap=cmap)
        axes[2,i].grid(False); axes[1,i].set_xticklabels([]); axes[1,i].set_yticklabels([]); plt.grid(False)

    # y軸ラベル
    axes[0,0].set_ylabel("Input image",fontsize=fontsize)
    axes[1,0].set_ylabel("True mask",fontsize=fontsize)
    axes[2,0].set_ylabel("Pred mask",fontsize=fontsize)
    

def Plot_perClassMetrics(
    history,
    metrix="dice_perclass",
    heatmap=True,
):
    if isinstance(history, tf.keras.callbacks.History):
        history = history.history
    num_classes = len(history[metrix][0])
    figsize= (7,4) if heatmap else (10,8)

    fig, axes = plt.subplots(2,1,squeeze=True,
                         figsize=figsize,
                         tight_layout=True,
                         facecolor="whitesmoke")
    
    if heatmap:
        axes[0].imshow(np.transpose(history[metrix]),vmin=0,vmax=1.0)
        axes[0].set_yticks(range(num_classes))        
        axes[1].imshow(np.transpose(history[f"val_{metrix}"]),vmin=0,vmax=1.0)
        axes[1].set_yticks(range(num_classes))      
    else:     
        axes[0].plot(history[metrix])
        axes[0].legend([f"class_{i}" for i in range(num_classes)],loc="right", bbox_to_anchor=(1.2, 0.5))
        axes[1].plot(history[f"val_{metrix}"])
        axes[1].legend([f"class_{i}" for i in range(num_classes)],loc="right", bbox_to_anchor=(1.2, 0.5))
    
    axes[0].set_title("train")
    axes[1].set_title("val")
    
    plt.show()


def CheckPrediction2(
    model,
    dataset,
    cmap="Spectral",
    num_classes=None,
    figsize=5,
    check_num:int=10
):
    if num_classes is not None:
        import matplotlib.colors as colors
        cm = plt.get_cmap(name=cmap, lut=num_classes-1) 
        cmap = colors.ListedColormap([(1.0, 1.0, 1.0)] + [ cm(i)[:3] for i in range(cm.N)])

    # datasetからcheck_num組の画像とラベルを取り出す
    dataset = dataset.unbatch().shuffle(buffer_size=1000)
    dataset = dataset.take(check_num)
    image_list = []; orig_image_list = []; truth_list = []
    for image, label in dataset:
        image_list.append(image)
        # one-hot対策
        if len(label.shape) > 2 and label.shape[-1] > 1:
            label = tf.argmax(label,axis=-1)
        truth_list.append(label)
        # 値の範囲を元に戻した明視野画像
        orig_image_list.append(tf.keras.preprocessing.image.array_to_img(image))

    # 予測用に画像リストを重ねてバッチ軸を作成
    images = tf.stack(image_list, axis=0)
    # 予測
    preds_list = model.predict(images)

    for i in range(len(preds_list)):
        preds = preds_list[i]
        # 最も確率が高いラベルへ
        preds = tf.argmax(preds, axis=-1)
    
        # 元画像/予測画像/正解画像 を並べて表示
        fig, axes = plt.subplots(ncols=check_num, nrows=3,
                                 figsize=(figsize*check_num,figsize*3),
                                 tight_layout=True)
        
        for i in range(check_num):
            axes[0,i].imshow(orig_image_list[i])
            axes[0,i].grid(False); axes[0,i].set_xticklabels([]); axes[0,i].set_yticklabels([]); plt.grid(False)
            axes[1,i].imshow(preds[i], vmin=0, vmax=num_classes, cmap=cmap)
            axes[1,i].grid(False); axes[1,i].set_xticklabels([]); axes[1,i].set_yticklabels([]); plt.grid(False)
            axes[2,i].imshow(truth_list[i], vmin=0, vmax=num_classes, cmap=cmap)
            axes[2,i].grid(False); axes[2,i].set_xticklabels([]); axes[2,i].set_yticklabels([]); plt.grid(False)        
    
        plt.show()