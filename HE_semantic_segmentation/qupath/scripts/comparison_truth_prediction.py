from pathlib import Path
import numpy as np
import pandas as pd
from PIL import Image
import tensorflow as tf

split_csv = Path("src/file_split.csv")
image_dir = Path("data/Train/Images")
label_dir = Path("data/Train/Labels")
model_path = Path("model_weights/model7_241111.h5")
output_dir = Path("qupath/val_prediction")
output_dir.mkdir(exist_ok=True)

seed = 1
n_tiles = 100

tissue_colors = {
    0:  "#000000",
    1:  "#a05aa0",
    2:  "#30d383",
    3:  "#83edaa",
    4:  "#ecffb3",
    5:  "#67da17",
    6:  "#44d2e5",
    7:  "#59220a",
    8:  "#f68d14",
    9:  "#808000",
    10: "#698c69",
    11: "#b299dc",
    12: "#ff00ff"
}

def hex_to_rgb(x):
    x = x.lstrip("#")
    return tuple(int(x[i:i+2], 16) for i in (0, 2, 4))

color_lut = np.zeros((256, 3), dtype=np.uint8)

for class_id, color in tissue_colors.items():
    color_lut[class_id] = hex_to_rgb(color)

df = pd.read_csv(split_csv)

val_df = df[df["split"] == "val"].copy()

print("Validation tiles:", len(val_df))

sample_df = val_df.sample(
    n=n_tiles,
    random_state=seed
).reset_index(drop=True)

print(sample_df)

sample_df.to_csv(
    output_dir / "selected_val_tiles.csv",
    index=False
)

model = tf.keras.models.load_model(
    model_path,
    compile=False
)

print("Model loaded")

for i, row in sample_df.iterrows():
    filename = row["filename"]
    image_path = image_dir / f"{filename}.jpg"
    label_path = label_dir / f"{filename}.png"
    print(f"{i+1}/100:", filename)
    image = np.array(Image.open(image_path).convert("RGB"), dtype=np.float32)
    truth = np.array(Image.open(label_path), dtype=np.uint8)
    x = image / 127.5 - 1.0
    x = np.expand_dims(x, axis=0)
    pred = model.predict(x, verbose=0)
    prediction = np.argmax(pred[0], axis=-1).astype(np.uint8)
    truth_rgb = color_lut[truth]
    prediction_rgb = color_lut[prediction]
    image_rgb = image.astype(np.uint8)
    gap = 10
    comparison = np.zeros((512, 512 * 3 + gap * 2, 3), dtype=np.uint8)
    comparison[:, 0:512] = image_rgb
    comparison[:, 512 + gap:512 + gap + 512] = truth_rgb
    comparison[:, 512 * 2 + gap * 2:] = prediction_rgb
    prefix = f"{i+1:02d}"
    Image.fromarray(image_rgb).save(output_dir / f"{prefix}_HE.png")
    Image.fromarray(truth_rgb).save(output_dir / f"{prefix}_truth.png")
    Image.fromarray(prediction_rgb).save(output_dir / f"{prefix}_prediction.png")
    Image.fromarray(prediction).save(output_dir / f"{prefix}_prediction_raw.png")
    Image.fromarray(comparison).save(output_dir / f"{prefix}_HE_truth_prediction.png")
    print("truth:", np.unique(truth), "prediction:", np.unique(prediction))

print("Done")
print("Output:", output_dir)
