import os
import csv
from pymol import cmd
from pymol.cgo import *


folder_path="/home/chiara/Desktop/COMENZE/Comenze_pipeline/enzymes_Review/pdb_files/pdb_aligned/"
output_csv="/home/chiara/Desktop/COMENZE/Comenze_pipeline/enzymes_Review/pdb_files/pdb_aligned/Structural_analysis/protein_dimension.csv"
def drawBox(obj_name, min_coord, max_coord):
    xmin, ymin, zmin = min_coord
    xmax, ymax, zmax = max_coord
    v=[[xmin, ymin, zmin], [xmax, ymin, zmin], [xmax, ymax, zmin], [xmin, ymax, zmin],
       [xmin, ymin, zmax], [xmax, ymin, zmax], [xmax, ymax, zmax], [xmin, ymax, zmax]]
    edges=[(0,1),(1,2),(2,3),(3,0),(4,5),(5,6),(6,7),(7,4),(0,4),(1,5),(2,6),(3,7)]

    box = [LINEWIDTH, 2.0, BEGIN, LINES, COLOR, 1.0, 0.0, 0.0]
    for e in edges:
        a, b = v[e[0]], v[e[1]]
        box += [VERTEX, a[0], a[1], a[2]]
        box += [VERTEX, b[0], b[1], b[2]]
    box += [END]

    cmd.load_cgo(box, "{}_bbox".format(obj_name))


def calculate_protein_dimension(obj_name):
    min_coord, max_coord = cmd.get_extent(obj_name)

    dx = max_coord [0] - min_coord [0]
    dy = max_coord [1] - min_coord [1]
    dz = max_coord [2] - min_coord [2]

    diagonal = (dx**2 + dy**2 + dz**2) ** 0.5

    print("{} dx: {} {} dy: {} {} dz: {} {}".format(
        obj_name,
        max_coord[0], min_coord[0],
        max_coord[1], min_coord[1],
        max_coord[2], min_coord[2]
    ))
    drawBox(obj_name, min_coord, max_coord)
    return dx, dy, dz, diagonal

def main():
    cmd.reinitialize()

    results = []

    for file in os.listdir(folder_path):
        if file.endswith(".pdb"):
            obj_name = os.path.splitext(file)[0]
            file_path = os.path.join(folder_path, file)

            cmd.load(file_path, obj_name)

            dx, dy, dz, diagonal = calculate_protein_dimension(obj_name)

            results.append((obj_name, dx, dy, dz, diagonal))
            cmd.hide("everything", obj_name)
            cmd.show("cartoon", obj_name)

    sorted_results = sorted(results, key=lambda x: x[4], reverse=True)

    with open(output_csv, "wb") as f:
        writer = csv.writer(f)
        writer.writerow(["protein", "dx", "dy", "dz", "diagonal"])
        writer.writerows(sorted_results)

    cmd.zoom("all")
    cmd.orient("all")
    

    print("\nNext 3 proteins after Top 5:")
    for r in sorted_results[5:8]:
        print("{} -> {:.2f} (dx={:.2f}, dy={:.2f}, dz={:.2f})".format(
            r[0], r[4], r[1], r[2], r[3]
        ))

    print("\nTop 3 smaller proteins:")

    for r in sorted_results[-3:]:

        print("{} -> {:.2f} (dx={:.2f}, dy={:.2f}, dz={:.2f})".format(
            r[0], r[4], r[1], r[2], r[3]
        ))
main()
