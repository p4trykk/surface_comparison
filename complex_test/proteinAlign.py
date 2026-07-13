
import os
import ast
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, PDBIO
import plotly.graph_objects as go
from IPython.display import display, HTML
import plotly.express as px
dir = 'pdb_from_pqr'
outDir = 'aligned_pdb'
csv = os.path.join(dir, 'tunnel_coords.csv')
os.makedirs(outDir, exist_ok=True)

targetPoint = np.array([-25.77, 1.56, -0.93])


def get_rotation_matrix(vec1, vec2):
    a = (vec1 / np.linalg.norm(vec1)).reshape(3)
    b = (vec2 / np.linalg.norm(vec2)).reshape(3)

    v = np.cross(a, b)
    c = np.dot(a, b)
    s = np.linalg.norm(v)

    if s < 1e-10:
        return np.eye(3)

    kmat = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])

    return np.eye(3) + kmat + kmat.dot(kmat) * ((1 - c) / (s ** 2))


df = pd.read_csv(csv)

parser = PDBParser(QUIET=True)
io = PDBIO()
aligned_results = []

maxPlotToShow = 2 #example

for index, row in df.iterrows():
    protein_name = row['Protein_name']

    tunnel_point = np.array([
        row['X'],
        row['Y'],
        row['Z']
    ], dtype=float)

    pdb_path = os.path.join(dir, f"{protein_name}.pdb")
    output_path = os.path.join(outDir, f"{protein_name}_aligned.pdb")

    if not os.path.exists(pdb_path):
        continue


    structure = parser.get_structure(protein_name, pdb_path)
    atoms = list(structure.get_atoms())
    coords = np.array([atom.get_coord() for atom in atoms])

    centroid = np.mean(coords, axis=0)

    centered_coords = coords - centroid
    centered_tunnel_point = tunnel_point - centroid

    target_axis = targetPoint  

    rot_matrix = get_rotation_matrix(centered_tunnel_point, target_axis)

    final_coords = np.dot(centered_coords, rot_matrix.T)
    final_tunnel_point = np.dot(centered_tunnel_point, rot_matrix.T)
    print(f"Final tunnel point: {final_tunnel_point}")

    final_vec_norm = final_tunnel_point / np.linalg.norm(final_tunnel_point)
    target_norm = target_axis / np.linalg.norm(target_axis)

    angle_after = np.degrees(np.arccos(np.clip(np.dot(final_vec_norm, target_norm), -1.0, 1.0)))

    print(f"angle after: {angle_after:.6f} deg")

    for atom, new_coord in zip(atoms, final_coords):
        atom.set_coord(new_coord)

    io.set_structure(structure)
    io.save(output_path)

    aligned_results.append({
        'name': protein_name,
        'coords': final_coords,
        'tunnel': final_tunnel_point
    })

print("done. ")
