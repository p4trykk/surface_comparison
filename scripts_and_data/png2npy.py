import os, re
import numpy as np
import matplotlib.pyplot as plt
from skimage.io import imread

INPUT_DIR  = "/content/drive/MyDrive/phD_interview/SOLENZ/dataToPCAplots/new_proteins_RADIUS_05" # do zmiany
OUTPUT_DIR = "/content/drive/MyDrive/phD_interview/SOLENZ/dataToPCAplots/new_proteins_RADIUS_05" # do zmiany x2
TARGET_SHAPE = None

os.makedirs(OUTPUT_DIR, exist_ok=True)

def parse_name(fname):
    stem = os.path.splitext(fname)[0]
    m = re.match(r"^(.+)_(rot\d{3})_bfactor_map_clean$", stem)
    return (m.group(1), m.group(2)) if m else (stem, None)

def png_to_potential_matrix(path):
    img = imread(path).astype(np.float32)
    if img.ndim == 2:
        img = img / 255.0 if img.max() > 1.0 else img
        return np.clip(img * 2.0 - 1.0, -1.0, 1.0).astype(np.float32)
    img = img[:, :, :3]
    if img.max() > 1.0: img = img / 255.0
    return np.clip(img[:,:,2] - img[:,:,0], -1.0, 1.0).astype(np.float32)


all_arrays = {}
files = sorted([f for f in os.listdir(INPUT_DIR) if f.lower().endswith('.png')])

for fname in files:
    path = os.path.join(INPUT_DIR, fname)
    protein, rot = parse_name(fname)
    mat = png_to_potential_matrix(path)
    ok_shape = (mat.shape == TARGET_SHAPE) if TARGET_SHAPE else True
    print(f'{fname}\n'
          f'  protein={protein} rot={rot} shape={mat.shape} '
          f'min/max=({mat.min():.3f},{mat.max():.3f}) ok={ok_shape}')
    base = os.path.splitext(fname)[0]
    np.save(os.path.join(OUTPUT_DIR, base + '.npy'), mat)
    all_arrays[base] = mat

combined_path = os.path.join(OUTPUT_DIR, 'all_maps.npz')
np.savez_compressed(combined_path, **all_arrays)
