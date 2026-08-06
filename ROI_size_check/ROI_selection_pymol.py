import csv
import re
from pathlib import Path

from pymol import cmd



CSV_PATH = Path(
    r"/home/chiara/Desktop/COMENZE/Comenze_pipeline/Surface_analysis/surfmap/output_SURFMAP_4OYY_hic_4rev_aligned_bfactor/smoothed_matrices/4OYY_hic_roi_residues.csv"
)

PDB_PATH = Path(
    r"/home/chiara/Desktop/COMENZE/Comenze_pipeline/enzymes_Review/pdb_files/pdb_aligned/4OYY_hic_4rev_aligned.pdb"
)

# Scegli "magenta" oppure "red"
RESIDUE_COLOR = "magenta"

# Nome della selezione che apparirà in PyMOL
SELECTION_NAME = "ROI_residues"

# True per mostrare anche gli sticks
SHOW_STICKS = False



# FUNZIONI


def create_object_name(pdb_path):
    """Crea un nome PyMOL valido dal nome del PDB."""
    name = re.sub(
        r"[^A-Za-z0-9_]",
        "_",
        pdb_path.stem,
    )

    if name[0].isdigit():
        name = "obj_" + name

    return name


def read_residues(csv_path):
    """
    Legge valori come:

    TYR_458_NA
    LYS_85_A

    Restituisce numero del residuo ed eventuale catena.
    """

    residues = []
    seen = set()

    with csv_path.open(
        "r",
        encoding="utf-8-sig",
        newline="",
    ) as file:

        reader = csv.DictReader(file)

        if "Residue" not in reader.fieldnames:
            raise ValueError(
                "Nel CSV non è presente la colonna 'Residue'."
            )

        for row in reader:
            label = row["Residue"].strip()

            if not label:
                continue

            parts = label.split("_")

            if len(parts) < 2:
                print(
                    f"Riga ignorata: formato non riconosciuto: {label}"
                )
                continue

            amino_acid = parts[0]
            residue_number = parts[1]

            chain = ""
            if len(parts) >= 3:
                chain = parts[2].strip()

            # NA significa che la catena non è specificata
            if chain.upper() in {
                "",
                "NA",
                "N/A",
                "NONE",
                "NAN",
            }:
                chain = ""

            # Evita residui duplicati
            key = (
                residue_number,
                chain,
            )

            if key not in seen:
                seen.add(key)

                residues.append(
                    (
                        amino_acid,
                        residue_number,
                        chain,
                    )
                )

    return residues



# ESECUZIONE


if not CSV_PATH.is_file():
    raise FileNotFoundError(
        f"CSV non trovato: {CSV_PATH}"
    )

if not PDB_PATH.is_file():
    raise FileNotFoundError(
        f"PDB non trovato: {PDB_PATH}"
    )

if RESIDUE_COLOR not in {
    "magenta",
    "red",
}:
    raise ValueError(
        "RESIDUE_COLOR deve essere 'magenta' oppure 'red'."
    )


# Carica il PDB
object_name = create_object_name(PDB_PATH)

if object_name in cmd.get_names("objects"):
    cmd.delete(object_name)

cmd.load(
    str(PDB_PATH),
    object_name,
)


# Legge i residui dal CSV
residues = read_residues(CSV_PATH)

selection_parts = []
found_residues = []
missing_residues = []


for amino_acid, residue_number, chain in residues:

    selection = (
        f"({object_name} "
        f"and polymer.protein "
        f"and resi {residue_number}"
    )

    if chain:
        selection += f" and chain {chain}"

    selection += ")"

    if cmd.count_atoms(selection) > 0:

        selection_parts.append(selection)

        found_residues.append(
            (
                amino_acid,
                residue_number,
                chain,
            )
        )

    else:

        missing_residues.append(
            (
                amino_acid,
                residue_number,
                chain,
            )
        )


if not selection_parts:
    raise ValueError(
        "Nessun residuo del CSV è stato trovato nel PDB."
    )


# Crea la selezione unica
cmd.select(
    SELECTION_NAME,
    " or ".join(selection_parts),
)


# Colora tutti gli atomi dei residui selezionati
cmd.color(
    RESIDUE_COLOR,
    SELECTION_NAME,
)


# Visualizzazione opzionale come sticks
if SHOW_STICKS:
    cmd.show(
        "sticks",
        SELECTION_NAME,
    )


# Centra la vista sulla selezione
cmd.zoom(
    SELECTION_NAME,
    buffer=5,
)

cmd.refresh()



# RISULTATO


print("\n" + "=" * 60)
print("SELEZIONE ROI COMPLETATA")
print("=" * 60)

print(
    f"Residui letti dal CSV: {len(residues)}"
)

print(
    f"Residui trovati nel PDB: {len(found_residues)}"
)

print(
    f"Atomi selezionati: {cmd.count_atoms(SELECTION_NAME)}"
)

print(
    f"Colore applicato: {RESIDUE_COLOR}"
)

print(
    f"Nome selezione PyMOL: {SELECTION_NAME}"
)


if missing_residues:

    print(
        f"\nResidui non trovati: {len(missing_residues)}"
    )

    for amino_acid, residue_number, chain in missing_residues:

        chain_text = (
            chain
            if chain
            else "non specificata"
        )

        print(
            f"- {amino_acid}_{residue_number} "
            f"catena {chain_text}"
        )

else:

    print(
        "\nTutti i residui sono stati trovati."
    )

print("=" * 60)