import argparse

parser = argparse.ArgumentParser()
parser.add_argument("model_path")
parser.add_argument("package_path")
parser.add_argument('image_dir_path')
parser.add_argument('--batch_size', default=50, type=int)
parser.add_argument('--tile_size', default=512, type=int)
parser.add_argument('--centering', default=True)
parser.add_argument('--save_dir', default="tmp_label")
parser.add_argument('--average_probability', default=False)
parser.add_argument('--label_splitsize', default=10000, type=int)
args = parser.parse_args()

import sys
sys.path.append(args.package_path)
from datasets import *
from prediction import *
from labelholder import *
from write_slide import *

def checkBool(arg):
    arg = arg.lower()
    if arg == "true":
        return True
    elif arg == "false":
        return False   

average_probability = checkBool(args.average_probability)
centering = checkBool(args.centering)

print("1. Read model")
model = tf.keras.models.load_model(args.model_path, compile=False)
print(model)

print(f"2. Create label holders from {args.image_dir_path}")
image_paths = glob(f"{args.image_dir_path}/*jpg")

labelholders = LabelHolders(image_dir=args.image_dir_path)
downsample=labelholders.labelholders[0].tile_info["d"]

print("3. Predict")
slide = labelholders.process(
    model, 
    batch_size=args.batch_size,
    centering=centering,
    average_probability=average_probability
    )

os.makedirs(args.save_dir, exist_ok=True)

min_y, max_y, min_x, max_x = getBoundingBox(slide)
slideArea = (max_y - min_y) * (max_x - min_x)
splitSize = args.label_splitsize
if splitSize > 0 and slideArea > splitSize*splitSize:
    saveSplitedLabel(
        slide, 
        tile_size = splitSize,
        downsample=downsample,
        slide_name=labelholders.slidename,
        save_dir=args.save_dir)
else:
    writeLabelSlide(
        slide,
        downsample=downsample,
        slide_name=labelholders.slidename,
        save_dir=args.save_dir
    )

print(f"4. label saved to {args.save_dir}")

del labelholders