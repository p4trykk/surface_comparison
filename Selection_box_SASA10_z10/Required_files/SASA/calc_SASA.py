import pandas as pd
from Bio.PDB import PDBParser
from Bio.PDB.SASA import ShrakeRupley
from itertools import chain
import os


parser = PDBParser(QUIET=True)
sr = ShrakeRupley()


# Values from Tien et al. 2013 - theoretical
MAX_SASA = {
    "A": 129,  # Alanine
    "R": 274,  # Arginine
    "N": 195,  # Asparagine
    "D": 193,  # Aspartate
    "C": 167,  # Cysteine
    "Q": 225,  # Glutamine
    "E": 223,  # Glutamate
    "G": 104,  # Glycine
    "H": 224,  # Histidine
    "I": 197,  # Isoleucine
    "L": 201,  # Leucine
    "K": 236,  # Lysine
    "M": 224,  # Methionine
    "F": 240,  # Phenylalanine
    "P": 159,  # Proline
    "S": 155,  # Serine
    "T": 172,  # Threonine
    "W": 285,  # Tryptophan
    "Y": 263,  # Tyrosine
    "V": 174   # Valine
}


def get_residue_code(residue):
    """
    Get the one-letter amino acid code for a residue.
    Args:
        residue: Residue object from Bio.PDB.
        
    Returns:
        str: One-letter code of the amino acid.
    """
    three_to_one = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
        'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V'
    }
    resname = residue.get_resname()
    return three_to_one.get(resname, None)

for filename in os.listdir("/home/chiara/Desktop/COMENZE/Comenze_pipeline/enzymes_Review/pdb_files/pdb_aligned/"):
    if filename.endswith(".pdb"):
        print(filename)
        protein = filename.split(".")[0]
        protein_sasa_dict = {'resid':[], 'resname':[], 'SASA':[], 'Relative SASA':[]}

        structure = parser.get_structure("protein", f"/home/chiara/Desktop/COMENZE/Comenze_pipeline/enzymes_Review/pdb_files/pdb_aligned/{filename}")
        # Compute SASA for the structure
        sr.compute(structure, level="R") 

        # Process the residues
        for model in structure:
            for chain in model:
                for residue in chain:
                    # Get SASA
                    sasa = residue.sasa
                    
                    # Get the one-letter amino acid code
                    resid = residue.id[1]
                    residue_code = get_residue_code(residue)

                    # Normalize by max SASA
                    max_sasa = MAX_SASA.get(residue_code)

                    if max_sasa is None:
                        continue

                    relative_sasa = (sasa / max_sasa) * 100

                    protein_sasa_dict['resid'].append(resid)
                    protein_sasa_dict['resname'].append(residue.resname)
                    protein_sasa_dict['SASA'].append(round(sasa, 2))
                    protein_sasa_dict['Relative SASA'].append(round(relative_sasa, 2))
        
        protein_sasa_df = pd.DataFrame(protein_sasa_dict)
        protein_sasa_df.to_csv(f'{protein}.csv', index=False)
        protein_sasa_df['Relative SASA'].to_csv(f'{protein}_RSA.txt', index=False, header=False)
    

