#!/usr/bin/env python3
"""
Batch conversion PQR -> PDB for the main enzyme dataset.

- Read coordinates/charges from the AMBER-generated PQR files.
- Write a standard fixed-column PDB.
- Store the PQR partial charge in the PDB B-factor field.
- Force chain A for every ATOM/HETATM record in this MAIN dataset.

"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List


# PATHS


PDB_PREPARATION_DIR = Path(
    "/home/chiara/Desktop/COMENZE/Comenze_pipeline/enzymes_Review/"
    "pdb_files/enzymes_4amber/output/preparation"
)

PQR_DIR = Path(
    "/home/chiara/Desktop/COMENZE/Comenze_pipeline/enzymes_Review/"
    "pdb_files/enzymes_4amber/output/enzymes_ready"
)

OUTPUT_DIR = Path(
    "/home/chiara/Desktop/COMENZE/Comenze_pipeline/enzymes_Review/"
    "pdb_files/enzymes_4amber/pqr_to_pdb_bfactor"
)

OUTPUT_CHAIN = "A"


# DATA

@dataclass
class PqrAtom:
    record: str
    serial: int
    name: str
    resname: str
    resseq: int
    x: float
    y: float
    z: float
    charge: float
    radius: float
    element: str


# PQR PARSING


def infer_element(atom_name: str, element_from_file: str = "") -> str:

    element = element_from_file.strip().upper()
    if element:
        return element[:2]

    name = atom_name.strip().upper()
    if not name:
        raise ValueError("Cannot infer element from an empty atom name.")

    # Protein hydrogen names can start with a digit or H.
    if name[0].isdigit():
        for ch in name:
            if ch.isalpha():
                return ch
    if name.startswith("H"):
        return "H"

    # Standard protein heavy atoms.
    if name.startswith("C"):
        return "C"
    if name.startswith("N"):
        return "N"
    if name.startswith("O"):
        return "O"
    if name.startswith("S"):
        return "S"
    if name.startswith("P"):
        return "P"

    # Conservative final fallback.
    letters = "".join(ch for ch in name if ch.isalpha())
    if not letters:
        raise ValueError(f"Cannot infer element for atom name {atom_name!r}")
    return letters[0]


def parse_pqr_atom(line: str) -> PqrAtom:

    record = line[0:6].strip()
    serial = int(line[6:11].strip())
    name = line[12:16].strip()
    resname = line[17:20].strip()


    if len(line) > 73:
        resseq = int(line[22:26].strip())
        x = float(line[30:38])
        y = float(line[38:46])
        z = float(line[46:54])
        charge = float(line[54:62])
        radius = float(line[62:70])
        raw_element = line[76:78].strip() if len(line) >= 78 else ""
    else:
        resseq = int(line[22:28].strip())
        x = float(line[28:38])
        y = float(line[38:46])
        z = float(line[46:54])
        charge = float(line[54:63])
        radius = float(line[63:72])
        raw_element = line[76:78].strip() if len(line) >= 78 else ""

    element = infer_element(name, raw_element)

    return PqrAtom(
        record=record,
        serial=serial,
        name=name,
        resname=resname,
        resseq=resseq,
        x=x,
        y=y,
        z=z,
        charge=charge,
        radius=radius,
        element=element,
    )



# PDB FORMATTING

def format_atom_name(atom_name: str, element: str) -> str:

    atom_name = atom_name.strip()[:4]
    element = element.strip().upper()

    if len(atom_name) == 4:
        return atom_name

    if atom_name and atom_name[0].isdigit():
        return f"{atom_name:<4s}"

    if len(element) == 1:
        return f" {atom_name:<3s}"[:4]

    return f"{atom_name:<4s}"[:4]


def format_pdb_atom(
    serial: int,
    atom: PqrAtom,
    chain: str = OUTPUT_CHAIN,
    occupancy: float = 1.00,
) -> str:
    """
    Column 22        : chain identifier
    Columns 61-66    : B-factor field -> PQR partial charge
    Columns 77-78    : element
    """
    element = infer_element(atom.name, atom.element)
    atom_name = format_atom_name(atom.name, element)

    line = (
        f"{atom.record:<6}"          #  1-6
        f"{serial:>5d}"              #  7-11
        f" "                         # 12
        f"{atom_name}"               # 13-16
        f" "                         # 17 altLoc
        f"{atom.resname:>3s}"        # 18-20
        f" "                         # 21
        f"{chain:1s}"                # 22 chain
        f"{atom.resseq:>4d}"         # 23-26
        f" "                         # 27 iCode
        f"   "                       # 28-30
        f"{atom.x:8.3f}"             # 31-38
        f"{atom.y:8.3f}"             # 39-46
        f"{atom.z:8.3f}"             # 47-54
        f"{occupancy:6.2f}"          # 55-60
        f"{atom.charge:6.2f}"        # 61-66 = B-factor / charge
        f"          "                # 67-76
        f"{element:>2s}"             # 77-78
        f"  "                        # 79-80
    )

    if len(line) != 80:
        raise RuntimeError(
            f"Internal PDB formatting error: expected 80 columns, got {len(line)}"
        )

    return line


def format_ter(serial: int, resname: str, resseq: int) -> str:
    """
    TER line with chain A. Kept simple because the main dataset is one chain.
    """
    return (
        f"TER   {serial:>5d}      "
        f"{resname:>3s} {OUTPUT_CHAIN:1s}{resseq:>4d}"
    )


# SOURCE / VALIDATION

def expected_preparation_pdb(pqr_path: Path) -> Path:
    """
        <stem>.pqr -> output/preparation/<stem>_final.pdb
    """
    return PDB_PREPARATION_DIR / f"{pqr_path.stem}_final.pdb"


def validate_output_pdb(path: Path) -> tuple[int, set[str]]:
    atom_count = 0
    chains = set()

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line_no, line in enumerate(handle, start=1):
            if not line.startswith(("ATOM  ", "HETATM")):
                continue

            atom_count += 1
            raw = line.rstrip("\n")

            if len(raw) != 80:
                raise RuntimeError(
                    f"{path.name}: line {line_no} has {len(raw)} columns, not 80."
                )

            chain = raw[21]
            chains.add(chain)

            if chain != OUTPUT_CHAIN:
                raise RuntimeError(
                    f"{path.name}: invalid chain {chain!r} at line {line_no}; "
                    f"expected only {OUTPUT_CHAIN!r}."
                )

            if not raw[76:78].strip():
                raise RuntimeError(
                    f"{path.name}: missing element in columns 77-78 at line {line_no}."
                )

    if atom_count == 0:
        raise RuntimeError(f"{path.name}: no ATOM/HETATM records were written.")

    return atom_count, chains


# CONVERSION

def convert(pdb_path: Path, pqr_path: Path, out_path: Path) -> None:
    """
    Final output is forced to A, selected chain A upstream.

    """
    if not pdb_path.exists():
        raise FileNotFoundError(
            f"Preparation PDB not found: {pdb_path}"
        )

    # Make sure the expected preparation PDB is not empty/corrupt.
    source_atom_count = sum(
        1
        for line in pdb_path.read_text(errors="replace").splitlines()
        if line.startswith("ATOM  ")
    )
    if source_atom_count == 0:
        raise RuntimeError(
            f"Preparation PDB contains no ATOM records: {pdb_path}"
        )

    pqr_lines = pqr_path.read_text(errors="replace").splitlines()

    pqr_atoms: List[PqrAtom] = [
        parse_pqr_atom(line)
        for line in pqr_lines
        if line.startswith(("ATOM  ", "HETATM"))
    ]

    if not pqr_atoms:
        raise RuntimeError(f"No ATOM/HETATM records in {pqr_path.name}")

    out_lines: List[str] = []

    # Re-number output atoms consecutively, as in the previous script.
    for new_serial, atom in enumerate(pqr_atoms, start=1):
        out_lines.append(
            format_pdb_atom(
                serial=new_serial,
                atom=atom,
                chain=OUTPUT_CHAIN,
            )
        )

    out_lines.append("END")
    out_path.write_text("\n".join(out_lines) + "\n", encoding="utf-8")

    validate_output_pdb(out_path)


# MAIN

def main() -> None:
    print("=" * 78)
    print("PQR -> PDB | charge in B-factor | FORCE CHAIN A")
    print("=" * 78)
    print(f"Preparation PDB dir : {PDB_PREPARATION_DIR}")
    print(f"PQR dir             : {PQR_DIR}")
    print(f"Output dir          : {OUTPUT_DIR}")
    print()

    if not PDB_PREPARATION_DIR.exists():
        raise FileNotFoundError(
            f"Preparation directory not found: {PDB_PREPARATION_DIR}"
        )

    if not PQR_DIR.exists():
        raise FileNotFoundError(f"PQR directory not found: {PQR_DIR}")

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    pqr_files = sorted(PQR_DIR.glob("*.pqr"))
    if not pqr_files:
        raise RuntimeError(f"No .pqr files found in {PQR_DIR}")

    converted = 0
    failed = 0
    excluded = 0

    for pqr_path in pqr_files:
        # 4CHA is a special A+B+C case and is intentionally handled separately.
        if "4CHA" in pqr_path.name.upper():
            print()
            print(f"SKIP SPECIAL CASE 4CHA: {pqr_path.name}")
            excluded += 1
            continue

        print()
        print("-" * 78)
        print(f"Processing: {pqr_path.name}")

        pdb_path = expected_preparation_pdb(pqr_path)
        out_path = OUTPUT_DIR / f"{pqr_path.stem}.pdb"

        print(f"  preparation PDB: {pdb_path.name}")
        print(f"  output PDB     : {out_path.name}")

        try:
            convert(pdb_path, pqr_path, out_path)

            atom_count, chains = validate_output_pdb(out_path)

            print(
                f"  DONE | atoms={atom_count} | chains={sorted(chains)} "
                f"| chain column 22 validated"
            )
            converted += 1

        except Exception as exc:
            print(f"  FAILED: {exc}")
            failed += 1

    print()
    print("=" * 78)
    print("FINISHED")
    print("=" * 78)
    print(f"Converted : {converted}")
    print(f"Failed    : {failed}")
    print(f"4CHA skip : {excluded}")
    print(f"Output    : {OUTPUT_DIR}")


if __name__ == "__main__":
    main()
