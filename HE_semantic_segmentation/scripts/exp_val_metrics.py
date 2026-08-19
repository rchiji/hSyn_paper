"""train / val / train+val の各タイル集合に model7 で推論し、クラス別 IoU / Dice / precision / recall と混同行列を算出する。

実行 (GPU): PATH に miniforge3/envs/tf210gpu/Library/bin を通した上で
    miniforge3/envs/tf210gpu/python.exe scripts/exp_val_metrics.py
出力: assets/{train,val,trainval}_confusion_matrix.csv, assets/{train,val,trainval}_class_metrics.csv
"""
import logging
import random
from glob import glob
from pathlib import Path

import numpy as np
from PIL import Image

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(message)s")
logger: logging.Logger = logging.getLogger(__name__)

BASE: Path = Path(__file__).resolve().parent.parent
DATA: Path = BASE / "data" / "Train"
ASSETS: Path = BASE / "assets"
MODEL_PATH: Path = BASE / "model_weights" / "model7_241111.h5"
N_GT: int = 13    # GT ラベルは 0–12
N_PRED: int = 14  # モデル出力は 14 クラス (class13 は GT に存在しない)
BATCH: int = 16

CLASS_NAMES: list[str] = [
    "Background", "Immune cells", "plasma", "Fibro(loose)",
    "Fibro(dense,regular)", "Fibro(dense,irregular)", "lining",
    "vessel", "vessel(large)", "adipose", "Stroma", "muscle", "RBC",
]


def evaluate_tiles(model, image_paths: list[str], mask_paths: list[str]) -> np.ndarray:
    """タイル集合を推論して混同行列 (GT 13 × pred 14) を返す。"""
    conf: np.ndarray = np.zeros((N_GT, N_PRED), dtype=np.int64)
    for b in range(0, len(image_paths), BATCH):
        imgs: np.ndarray = np.stack([
            np.asarray(Image.open(p), dtype=np.float32) / 127.5 - 1.0
            for p in image_paths[b:b + BATCH]
        ])  # (B,512,512,3)
        preds: np.ndarray = model.predict(imgs, verbose=0).argmax(-1)  # (B,512,512)
        gts: np.ndarray = np.stack([
            np.asarray(Image.open(p)) for p in mask_paths[b:b + BATCH]
        ])  # (B,512,512)
        conf += np.bincount(
            (gts.ravel().astype(np.int64) * N_PRED + preds.ravel()),
            minlength=N_GT * N_PRED,
        ).reshape(N_GT, N_PRED)
        if (b // BATCH) % 20 == 0:
            logger.info("batch %d / %d", b // BATCH, len(image_paths) // BATCH)
    return conf


def save_metrics(conf: np.ndarray, split: str) -> None:
    """混同行列からクラス別メトリクスを計算して CSV 保存・ログ出力する。"""
    np.savetxt(ASSETS / f"{split}_confusion_matrix.csv", conf, fmt="%d", delimiter=",",
               header=",".join(CLASS_NAMES + ["unused13"]), comments="")
    tp: np.ndarray = np.diag(conf[:, :N_GT]).astype(np.float64)  # (13,)
    fp: np.ndarray = conf.sum(0)[:N_GT] - tp
    fn: np.ndarray = conf.sum(1) - tp
    iou: np.ndarray = tp / np.maximum(tp + fp + fn, 1)
    dice: np.ndarray = 2 * tp / np.maximum(2 * tp + fp + fn, 1)
    precision: np.ndarray = tp / np.maximum(tp + fp, 1)
    recall: np.ndarray = tp / np.maximum(tp + fn, 1)
    support: np.ndarray = conf.sum(1)

    with open(ASSETS / f"{split}_class_metrics.csv", "w", encoding="utf-8") as f:
        f.write("class,name,IoU,Dice,precision,recall,support_px\n")
        for k in range(N_GT):
            f.write(f'{k},"{CLASS_NAMES[k]}",{iou[k]:.4f},{dice[k]:.4f},'
                    f"{precision[k]:.4f},{recall[k]:.4f},{int(support[k])}\n")

    logger.info("[%s] mean IoU (all): %.4f / (excl. background): %.4f", split, iou.mean(), iou[1:].mean())
    logger.info("[%s] pixel accuracy: %.4f", split, tp.sum() / conf.sum())
    for k in range(N_GT):
        logger.info("[%s] %-24s IoU=%.3f Dice=%.3f P=%.3f R=%.3f", split, CLASS_NAMES[k], iou[k], dice[k], precision[k], recall[k])


def main() -> None:
    import tensorflow as tf

    ipaths: list[str] = sorted(glob(str(DATA / "Images" / "*jpg")))
    mpaths: list[str] = sorted(glob(str(DATA / "Labels" / "*png")))
    assert len(ipaths) == len(mpaths) == 2690

    # PrepareTrainValPath (seed=1, train_ratio=0.8) と同一ロジックで分割を再現
    random.seed(1)
    train_idx: set[int] = set(random.sample(range(len(ipaths)), int(len(ipaths) * 0.8)))
    splits: dict[str, tuple[list[str], list[str]]] = {
        "train": ([p for i, p in enumerate(ipaths) if i in train_idx],
                  [p for i, p in enumerate(mpaths) if i in train_idx]),
        "val": ([p for i, p in enumerate(ipaths) if i not in train_idx],
                [p for i, p in enumerate(mpaths) if i not in train_idx]),
    }

    model = tf.keras.models.load_model(MODEL_PATH, compile=False)
    logger.info("model loaded: input %s output %s", model.input_shape, model.output_shape)

    confs: dict[str, np.ndarray] = {}
    for split, (ips, mps) in splits.items():
        logger.info("[%s] %d tiles", split, len(ips))
        confs[split] = evaluate_tiles(model, ips, mps)
        save_metrics(confs[split], split)

    # train+val は両混同行列の和 (同一タイル集合の合併なので再推論不要)
    save_metrics(confs["train"] + confs["val"], "trainval")


if __name__ == "__main__":
    main()
