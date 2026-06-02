import os
import subprocess

INPUT_DIR = "."
LOG_DIR = "surfmap_logs"

os.makedirs(LOG_DIR, exist_ok=True)

print("Start >>>\n")

for file in os.listdir(INPUT_DIR):
    if not file.endswith(".pdb"):
        continue

    pdb_path = os.path.join(INPUT_DIR, file)
    base = file.replace(".pdb", "")

    log_file = os.path.join(LOG_DIR, f"{base}.log")

    print(f"{file}")

    cmd = [
        "surfmap",
        "-pdb", pdb_path,
        "-tomap", "bfactor",
        "--png",
        "--keep",
        "--clean-map",
        "-s", "2",
        "--bfactor-normalize"
    ]

    try:
        with open(log_file, "w") as log:
            result = subprocess.run(cmd, stdout=log, stderr=log)

        if result.returncode != 0:
            print(f" {log_file}\n")
            continue

    except Exception as e:
        print(f"{e}")
        continue

print("FINISHED :| ")
