from pymol import cmd
from pymol.cgo import *
import os
import pandas as pd
import numpy as np


# PATHS
path_SASA = "/home/chiara/Desktop/COMENZE/Comenze_pipeline/enzymes_Review/SASA"
path_pdb = "/home/chiara/Desktop/COMENZE/Comenze_pipeline/enzymes_Review/pdb_files/pdb_aligned"
path_coords = "/home/chiara/Desktop/COMENZE/Comenze_pipeline/enzymes_Review/pdb_files/pdb_aligned/Structural_analysis/aligned_proteins_tunnel_coords.txt"
path_dims = "/home/chiara/Desktop/COMENZE/Comenze_pipeline/enzymes_Review/pdb_files/pdb_aligned/Structural_analysis/protein_dimension.csv"
output_csv_path = "/home/chiara/Desktop/COMENZE/Comenze_pipeline/enzymes_Review/pdb_files/pdb_aligned/Structural_analysis"


# PARAMETERS
sasa_threshold = 10.0
box_front = 10.0
box_back = 5.0  # Å, depth of the selection box
tunnel_vec = np.array([
    -25.77,
     1.56,
    -0.93
])

tunnel_vec = tunnel_vec / np.linalg.norm(tunnel_vec)

# tunnel coordinates
coords_list = []

with open(path_coords, "r") as f:
    for line in f:

        line = line.strip()

        if not line:
            continue

        protein, coords = line.split(",[")

        coords = coords.replace("]", "")

        x, y, z = coords.split(",")

        coords_list.append({
            "protein": protein.strip(),
            "x": float(x),
            "y": float(y),
            "z": float(z)
        })

coords_df = pd.DataFrame(coords_list)

# protein dimensions
dim_df = pd.read_csv(path_dims)
dim_df.columns = dim_df.columns.str.strip()
dim_df["protein"] = dim_df["protein"].astype(str).str.strip()

all_selected_residues = []


# LOOP PDB FILES
for pdb_file in os.listdir(path_pdb):

    if not pdb_file.endswith(".pdb"):
        continue

    protein = pdb_file.replace(".pdb", "")
    protein_base = protein.replace("_aligned", "")

    print(f"\nProcessing {protein}")

    cmd.load(os.path.join(path_pdb, pdb_file), protein)

    # SASA FILE

    sasa_file = os.path.join(path_SASA, protein + ".csv")

    if not os.path.exists(sasa_file):
        print(f"Missing SASA file for {protein}")
        continue

    sasa_df = pd.read_csv(sasa_file)
    sasa_df.columns = sasa_df.columns.str.strip()

    sasa_df_filtered = sasa_df[
        sasa_df["Relative SASA"] > sasa_threshold
    ]

    # TUNNEL COORDINATES

    coord_row = coords_df[
        coords_df["protein"] == protein_base
    ]

    if coord_row.empty:
        print(f"No coordinates for {protein}")
        continue

    x = float(coord_row["x"].values[0])
    y = float(coord_row["y"].values[0])
    z = float(coord_row["z"].values[0])

    center_front = np.array([x, y, z])

    v = tunnel_vec

    tmp = np.array([0, 0, 1])

    if abs(np.dot(v, tmp)) > 0.95:
        tmp = np.array([0, 1, 0])

    u = np.cross(v, tmp)
    u = u / np.linalg.norm(u)

    w = np.cross(v, u)
    w = w / np.linalg.norm(w)

    # PROTEIN DIMENSIONS

    dim_row = dim_df[
        dim_df["protein"] == protein
    ]

    if dim_row.empty:
        print(f"No dimension for {protein}")
        continue

    dx = float(dim_row["dx"].values[0])
    dy = float(dim_row["dy"].values[0])

    print(
        f"Box dimensions = "
        f"{dx:.2f} x {dy:.2f} x {(box_front + box_back):.2f}"
    )
 
    # BOX GEOMETRY

    # Tunnel point = center of front face
    # Box extends toward negative Z

    half_dx = dx / 2.0
    half_dy = dy / 2.0

    front_center = np.array([x, y, z])

    front_face_center = front_center + v * box_front
    back_face_center = front_center - v * box_back
    # VISUALIZE BOX

    box_name = f"{protein}_selection_box"

        # FRONT FACE

    f1 = front_face_center - half_dx * u - half_dy * w
    f2 = front_face_center + half_dx * u - half_dy * w
    f3 = front_face_center + half_dx * u + half_dy * w
    f4 = front_face_center - half_dx * u + half_dy * w

        # BACK FACE

    b1 = back_face_center - half_dx * u - half_dy * w
    b2 = back_face_center + half_dx * u - half_dy * w
    b3 = back_face_center + half_dx * u + half_dy * w
    b4 = back_face_center - half_dx * u + half_dy * w

    box_obj = [

        LINEWIDTH, 4.0,

        ALPHA, 1.0,

        COLOR, 0.0, 0.0, 1.0,

        BEGIN, LINES,
    # front

        VERTEX, *f1,
        VERTEX, *f2,

        VERTEX, *f2,
        VERTEX, *f3,

        VERTEX, *f3,
        VERTEX, *f4,

        VERTEX, *f4,
        VERTEX, *f1,

    # back

        VERTEX, *b1,
        VERTEX, *b2,

        VERTEX, *b2,
        VERTEX, *b3,

        VERTEX, *b3,
        VERTEX, *b4,

        VERTEX, *b4,
        VERTEX, *b1,

    # connections

        VERTEX, *f1,
        VERTEX, *b1,

        VERTEX, *f2,
        VERTEX, *b2,

        VERTEX, *f3,
        VERTEX, *b3,

        VERTEX, *f4,
        VERTEX, *b4,

        END
    ]

    cmd.load_cgo(box_obj, box_name)

 
    # TUNNEL POINT

    tunnel_point = f"{protein}_tunnel_point"

    cmd.pseudoatom(
        tunnel_point,
        pos=[x, y, z]
    )

    direction_point = f"{protein}_direction"
    direction_coord = (front_center - v * box_back).tolist()

    cmd.pseudoatom(
        direction_point,
        pos=direction_coord
    )
    cmd.show("spheres", direction_point)
    cmd.color("green", direction_point)

    cmd.show("spheres", tunnel_point)
    cmd.color("yellow", tunnel_point)

 
    # SELECTIONS

    selection_box = f"{protein}_BOX"
    selection_sasa = f"{protein}_BOX_SASA10"

    cmd.select(selection_box, "none")
    cmd.select(selection_sasa, "none")

    selected_residues_all = []
    selected_residues_sasa = []

    model = cmd.get_model(
        f"{protein} and name CA"
    )

    for atom in model.atom:

        resid = int(atom.resi)
        resn = atom.resn

        coord = np.array(atom.coord)

        r = coord - front_center

        local_u = np.dot(r, u)
        local_w = np.dot(r, w)
        local_v = np.dot(r, -v)

        if resid == int(model.atom[0].resi):
            print(
                f"Resid {resid}: "
                f"u={local_u:.2f} "
                f"w={local_w:.2f} "
                f"v={local_v:.2f}"
            )

        inside_box = (

            abs(local_u) <= half_dx and
            abs(local_w) <= half_dy and
            -box_front <= local_v <= box_back

        )

        if not inside_box:
            continue

        residue_selection = (
            f"({protein} and resi {resid})"
        )

        cmd.select(
            selection_box,
            f"{selection_box} or {residue_selection}"
        )

        selected_residues_all.append(
            f"{resn}{resid}"
        )

        sasa_row = sasa_df_filtered[
            sasa_df_filtered["resid"] == resid
        ]

        if sasa_row.empty:
            continue

        sasa_value = float(
            sasa_row["Relative SASA"].values[0]
        )

        cmd.select(
            selection_sasa,
            f"{selection_sasa} or {residue_selection}"
        )

        selected_residues_sasa.append(
            f"{resn}{resid}({sasa_value:.1f})"
        )


    # VISUALIZATION

    cmd.show("sticks", selection_box)
    cmd.color("orange", selection_box)

    cmd.show("sticks", selection_sasa)
    cmd.color("red", selection_sasa)


    # GROUP EVERYTHING

    group_name = f"{protein}_GROUP"

    cmd.group(
        group_name,
        " ".join([
            box_name,
            tunnel_point,
            direction_point,
            selection_box,
            selection_sasa
        ])
    )


    # SAVE RESULTS

    all_selected_residues.append({

        "protein": protein,

        "all_residues_in_box":
            "; ".join(selected_residues_all),

        "sasa_residues_in_box":
            "; ".join(selected_residues_sasa)

    })

    final_df = pd.DataFrame(
        all_selected_residues
    )

    final_output = os.path.join(
        output_csv_path,
        "Selection_SASA10_z10_Box.csv"
    )

    final_df.to_csv(
        final_output,
        index=False
    )

    print(f"{protein} completed")

print("\nDONE")
print("Saved final output:")
print(final_output)