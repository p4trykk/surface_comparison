#usage:: python3 pqr_to_pdb_bfactor.py input.pdb input.pqr output.pdb
from __future__ import annotations
import argparse
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple


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


@dataclass
class PdbAtomInfo:
    chain: str
    resname: str
    resseq: int
    i_code: str


def parse_pqr_atom(line: str) -> PqrAtom:
    record = line[0:6].strip()
    serial = int(line[6:11])
    name = line[12:16].strip()
    resname = line[17:20].strip()
    resseq = int(line[22:26])
    x = float(line[30:38])
    y = float(line[38:46])
    z = float(line[46:54])
    charge = float(line[54:62])
    radius = float(line[62:70])
    element = line[76:78].strip() if len(line) >= 78 else ""
    if not element:
        element = "".join(ch for ch in name if ch.isalpha())[:2].strip().upper()
    return PqrAtom(record, serial, name, resname, resseq, x, y, z, charge, radius, element)


def parse_pdb_metadata(pdb_lines: List[str]) -> Tuple[Dict[int, PdbAtomInfo], Set[int]]:
    resinfo: Dict[int, PdbAtomInfo] = {}
    ter_after: Set[int] = set()

    last_atom_resseq: Optional[int] = None
    for i, line in enumerate(pdb_lines):
        rec = line[0:6].strip()
        if rec in ("ATOM", "HETATM"):
            resseq = int(line[22:26])
            chain = line[21].strip() or " "
            resname = line[17:20].strip()
            i_code = line[26].strip() or " "
            resinfo[resseq] = PdbAtomInfo(chain=chain, resname=resname, resseq=resseq, i_code=i_code)
            last_atom_resseq = resseq
        elif rec == "TER" and last_atom_resseq is not None:
            ter_after.add(last_atom_resseq)

    return resinfo, ter_after


def format_pdb_atom(serial: int, atom: PqrAtom, chain: str, occupancy: float = 1.00) -> str:
    name = atom.name[:4]
    return (
        f"{atom.record:<6}{serial:>5} "
        f"{name:>4s}"
        f" "
        f"{atom.resname:>3s} "
        f"{chain:1s}"
        f"{atom.resseq:>4d}"
        f"    "
        f"{atom.x:8.3f}{atom.y:8.3f}{atom.z:8.3f}"
        f"{occupancy:6.2f}{atom.charge:6.2f}"
        f"          "
        f"{atom.element:>2s}"
    )


def format_ter(serial: int, resname: str, chain: str, resseq: int) -> str:
    return f"TER   {serial:>5d}      {resname:>3s} {chain:1s}{resseq:>4d}"


def convert(pdb_path: Path, pqr_path: Path, out_path: Path) -> None:
    pdb_lines = pdb_path.read_text(errors="replace").splitlines()
    pqr_lines = pqr_path.read_text(errors="replace").splitlines()

    resinfo, ter_after = parse_pdb_metadata(pdb_lines)

    pqr_atoms: List[PqrAtom] = [parse_pqr_atom(l) for l in pqr_lines if l.startswith(("ATOM", "HETATM"))]
    if not pqr_atoms:
        raise ValueError("No ATOM/HETATM in pqr")

    out_lines: List[str] = []
    serial = 1
    prev_resseq: Optional[int] = None
    prev_resname: Optional[str] = None
    prev_chain: Optional[str] = None

    for atom in pqr_atoms:
        info = resinfo.get(atom.resseq)
        chain = info.chain if info else " "
        if prev_resseq is not None and atom.resseq != prev_resseq:
            if prev_resseq in ter_after and prev_resname is not None and prev_chain is not None:
                out_lines.append(format_ter(serial, prev_resname, prev_chain, prev_resseq))
                serial += 1

        out_lines.append(format_pdb_atom(serial, atom, chain))
        serial += 1
        prev_resseq = atom.resseq
        prev_resname = atom.resname
        prev_chain = chain

    if prev_resseq is not None and prev_resseq in ter_after and prev_resname is not None and prev_chain is not None:
        out_lines.append(format_ter(serial, prev_resname, prev_chain, prev_resseq))
        serial += 1

    out_lines.append("END")
    out_path.write_text("\n".join(out_lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("input_pdb", type=Path)
    parser.add_argument("input_pqr", type=Path)
    parser.add_argument("output_pdb", type=Path)
    args = parser.parse_args()
    convert(args.input_pdb, args.input_pqr, args.output_pdb)

if __name__ == "__main__":
    main()
