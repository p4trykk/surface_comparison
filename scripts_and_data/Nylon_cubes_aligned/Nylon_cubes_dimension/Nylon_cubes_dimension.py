import os
import csv
from pymol import cmd
from pymol.cgo import *

# ============================================================
# PATHS
# ============================================================

folder_path = "/home/chiara/Desktop/COMENZE/Surfmatch_pipeline/Nylon_surface/Nylon_cube/final_structure"

# CSV salvato nella stessa cartella dei PDB
output_csv = os.path.join(folder_path, "nylon_cubes_dimension.csv")


# ============================================================
# DRAW BOUNDING BOX
# ============================================================

def drawBox(obj_name, min_coord, max_coord):
    xmin, ymin, zmin = min_coord
    xmax, ymax, zmax = max_coord

    vertices = [
        [xmin, ymin, zmin],
        [xmax, ymin, zmin],
        [xmax, ymax, zmin],
        [xmin, ymax, zmin],
        [xmin, ymin, zmax],
        [xmax, ymin, zmax],
        [xmax, ymax, zmax],
        [xmin, ymax, zmax]
    ]

    edges = [
        (0, 1), (1, 2), (2, 3), (3, 0),
        (4, 5), (5, 6), (6, 7), (7, 4),
        (0, 4), (1, 5), (2, 6), (3, 7)
    ]

    box = [LINEWIDTH, 2.0, BEGIN, LINES, COLOR, 1.0, 0.0, 0.0]

    for e in edges:
        a, b = vertices[e[0]], vertices[e[1]]
        box += [VERTEX, a[0], a[1], a[2]]
        box += [VERTEX, b[0], b[1], b[2]]

    box += [END]

    cmd.load_cgo(box, "{}_bbox".format(obj_name))


# ============================================================
# CALCULATE NYLON CUBE DIMENSIONS
# ============================================================

def calculate_cube_dimension(obj_name):
    min_coord, max_coord = cmd.get_extent(obj_name)

    dx = max_coord[0] - min_coord[0]
    dy = max_coord[1] - min_coord[1]
    dz = max_coord[2] - min_coord[2]

    diagonal = (dx**2 + dy**2 + dz**2) ** 0.5

    print(
        "{} -> dx={:.3f} A, dy={:.3f} A, dz={:.3f} A, diagonal={:.3f} A".format(
            obj_name, dx, dy, dz, diagonal
        )
    )

    drawBox(obj_name, min_coord, max_coord)

    return dx, dy, dz, diagonal


# ============================================================
# SAVE CSV
# ============================================================

def save_csv(results):
    """
    Salva progressivamente tutti i risultati disponibili.
    I risultati vengono ordinati dalla diagonale maggiore alla minore.
    """

    sorted_results = sorted(results, key=lambda x: x[4], reverse=True)

    with open(output_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)

        writer.writerow([
            "nylon_cube",
            "dx_A",
            "dy_A",
            "dz_A",
            "diagonal_A"
        ])

        writer.writerows(sorted_results)

    return sorted_results


# ============================================================
# MAIN
# ============================================================

def main():
    cmd.reinitialize()

    if not os.path.isdir(folder_path):
        raise FileNotFoundError(
            "La cartella non esiste:\n{}".format(folder_path)
        )

    # Legge SOLO i PDB che terminano con _clean.pdb
    pdb_files = sorted(
        [
            f for f in os.listdir(folder_path)
            if f.lower().endswith("_clean.pdb")
        ]
    )

    if not pdb_files:
        raise RuntimeError(
            "Non ho trovato file che terminano con _clean.pdb nella cartella:\n{}".format(
                folder_path
            )
        )

    results = []

    # Crea immediatamente il CSV con l'intestazione
    save_csv(results)

    print("\n============================================================")
    print("NYLON CUBE DIMENSION ANALYSIS")
    print("============================================================")
    print("PDB _clean trovati: {}".format(len(pdb_files)))
    print("Cartella PDB:")
    print(folder_path)
    print("CSV di output:")
    print(output_csv)
    print("============================================================\n")

    for i, file in enumerate(pdb_files, start=1):

        obj_name = os.path.splitext(file)[0]
        file_path = os.path.join(folder_path, file)

        print("[{}/{}] {}".format(i, len(pdb_files), obj_name))

        try:
            cmd.load(file_path, obj_name)

            dx, dy, dz, diagonal = calculate_cube_dimension(obj_name)

            results.append(
                (obj_name, dx, dy, dz, diagonal)
            )

            # Salva subito dopo ogni struttura
            save_csv(results)

            cmd.hide("everything", obj_name)
            cmd.show("sticks", obj_name)

        except Exception as e:
            print("ERRORE con {}: {}".format(file, e))
            print("Continuo con il prossimo PDB...\n")
            continue

    sorted_results = save_csv(results)

    cmd.zoom("all")
    cmd.orient("all")

    print("\n============================================================")
    print("ANALISI COMPLETATA")
    print("Cubi di nylon elaborati correttamente: {}".format(len(sorted_results)))
    print("CSV salvato in:")
    print(output_csv)
    print("============================================================")

    if sorted_results:
        print("\nDimensioni calcolate:")

        for r in sorted_results:
            print(
                "{} -> dx={:.2f} A, dy={:.2f} A, dz={:.2f} A, diagonal={:.2f} A".format(
                    r[0], r[1], r[2], r[3], r[4]
                )
            )


main()
