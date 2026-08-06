import pandas as pd
from collections import defaultdict
import math

TXT_FILE = "1CEX_cutfv_4rev_aligned_bfactor_smoothed_matrix.txt"
CENTER_X = 90
CENTER_Y = 45
RADIUS = 25
OUTPUT_CSV = "1CEX_roi_residues.csv"

# df = pd.read_csv(
#     TXT_FILE,
#     sep=r"\s+",
#     header=None,
#     names=["txt_x", "txt_y", "value", "residues"]
# )
records = []

with open(TXT_FILE) as f:

    for line in f:

        line = line.strip()

        if not line:
            continue

        parts = line.split(maxsplit=3)

        if len(parts) < 4:
            continue

        if parts[0].lower() == "absc":
            continue
        if not parts[0].isdigit():
            continue

        txt_x = int(parts[0])
        txt_y = int(parts[1])

        value = parts[2]
        residues = parts[3]

        records.append([txt_x, txt_y, value, residues])

df = pd.DataFrame(records, columns=["txt_x","txt_y","value","residues"])

df["x"] = df["txt_x"] // 2
df["y"] = df["txt_y"] // 2

radius2 = RADIUS ** 2

roi = df[
    (df["x"] - CENTER_X) ** 2 +
    (df["y"] - CENTER_Y) ** 2 <= radius2
].copy()

print(f"Liczba pikseli w ROI: {len(roi)}")

stats = defaultdict(lambda: {
    "pixels": 0,
    "values": []
})

for _, row in roi.iterrows():

    residues = str(row["residues"]).strip()

    if residues == "NA":
        continue

    value = row["value"]

    try:
        value = float(value)
    except:
        value = float("nan")

    for residue in residues.split(","):

        residue = residue.strip()

        if residue == "":
            continue

        stats[residue]["pixels"] += 1

        if not math.isnan(value):
            stats[residue]["values"].append(value)

results = []

for residue, data in stats.items():

    values = data["values"]

    if len(values) > 0:
        mean_value = sum(values) / len(values)
        min_value = min(values)
        max_value = max(values)
    else:
        mean_value = float("nan")
        min_value = float("nan")
        max_value = float("nan")

    results.append({"Residue": residue, "Pixel_count": data["pixels"], "Mean_value": mean_value, "Min_value": min_value, "Max_value": max_value})

results = pd.DataFrame(results)

results = results.sort_values(by="Pixel_count", ascending=False)

results.to_csv(OUTPUT_CSV, index=False)

print(results)

print()
print(f"found: {len(results)} aminoacids")
