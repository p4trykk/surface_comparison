#!/usr/bin/env python3

"""
Standardized batch CAVER pipeline using the terminal catalytic hydrogen and its XYZ coordinates.

INPUT
-----
A CSV containing at least:
    protein_name
    pdb_code

Starting point:
    all proteins -> chain A
    4CHA         -> chains A+B+C, catalytic starting chain C

If the same protein name appears more than once in the CSV, the PDB/source
code is appended to keep folders unique, e.g.:
    CAVER_Terminal/UMGSP2_9FZW/
    CAVER_Terminal/UMGSP2_8WDW/

Prints all failed proteins with the reason at the end.
"""

from __future__ import annotations

import csv
import re
import shutil
import subprocess
from collections import Counter
from pathlib import Path



# SETTINGS

# CSV: put it in the same directory as this script, or change this path.
CSV_FILE = Path(__file__).resolve().with_name("caver_starting_points.csv")

# Folder containing the PDB files converted from PQR, all ending in _4rev.pdb.
PDB_DIR = Path(
    "/home/chiara/Desktop/COMENZE/Comenze_pipeline/enzymes_Review/pdb_files/enzymes_4amber/pqr_to_pdb_bfactor"
)

# All standardized terminal CAVER results go here.
CAVER_TERMINAL_ROOT = Path(
    "/home/chiara/Desktop/COMENZE/Comenze_pipeline/enzymes_Review/Caver_clean/CAVER_Terminal"
)

# CAVER installation already tested from terminal.
CAVER_HOME = Path("/home/chiara/software/caver_3.0/caver")

JAVA_XMX = "4000m"

# Parameters 
PROBE_RADIUS = 1.4
SHELL_RADIUS = 10.0
SHELL_DEPTH = 4.0
FRAME_WEIGHTING_COEFFICIENT = 1.0
FRAME_CLUSTERING_THRESHOLD = 1.0

# Maximum distance (Angstrom) allowed by CAVER from starting point
MAX_DISTANCE = 3.0

DESIRED_RADIUS = 5.0


OVERWRITE = True

# For a first test, for example:
# ONLY_PROTEINS = {"ReU"}
# For the complete dataset:
ONLY_PROTEINS = None



# RESIDUE NAMES

CANONICAL_AA = {
    "ALA", "ARG", "ASN", "ASP", "CYS",
    "GLN", "GLU", "GLY", "HIS", "ILE",
    "LEU", "LYS", "MET", "PHE", "PRO",
    "SER", "THR", "TRP", "TYR", "VAL",
}

# CAVER doen't recognize Amber residues
AMBER_TO_STANDARD = {
    "HID": "HIS",
    "HIE": "HIS",
    "HIP": "HIS",
    "CYX": "CYS",
    "CYM": "CYS",
    "ASH": "ASP",
    "GLH": "GLU",
    "LYN": "LYS",
}

# Exact terminal hydrogen used as CAVER starting point.
TERMINAL_CATALYTIC_HYDROGEN = {
    "SER": "HG",
    "THR": "HG1",
    "CYS": "HG",
}

# Manual starting-point coordinates for specific buried catalytic sites.
# ONLY these proteins use the coordinates below.
MANUAL_STARTING_POINTS = {
    "NylB": (
        50.8446044921875,
        9.462358474731445,
        15.43001937866211,
    ),
    "UMGSP3": (
        -10.204702377319336,
        -40.20881652832031,
        18.575401306152344,
    ),
    "NfPolyA": (
        -3.6654775142669678,
        -11.782846450805664,
        6.403104782104492,
    ),
}


# CHAIN SELECTION

# Default for the whole dataset: chain A only.
DEFAULT_ALLOWED_CHAINS = ("A",)
DEFAULT_STARTING_CHAIN = "A"

# Special case:
# 4CHA (Chymotrypsin) must retain chains A+B+C.
# The catalytic serine is searched specifically in chain C.
PDB_CHAIN_OVERRIDES = {
    "4CHA": ("A", "B", "C"),
}

PDB_STARTING_CHAIN_OVERRIDES = {
    "4CHA": "C",
}


def parse_chain_list(value: str) -> tuple[str, ...]:
    """
    Parse optional CSV chain strings such as:
        A
        A+B+C
        A,B,C
        A B C
    """
    value = (value or "").strip()
    if not value:
        return tuple()

    chains = [
        item.strip()
        for item in re.split(r"[+,;\s]+", value)
        if item.strip()
    ]

    return tuple(dict.fromkeys(chains))


def get_allowed_chains(row: dict) -> tuple[str, ...]:
    """
    Priority:
      1. optional CSV column 'chains'
      2. PDB-specific override (4CHA -> A+B+C)
      3. default chain A
    """
    csv_chains = parse_chain_list(row.get("chains", "") or "")
    if csv_chains:
        return csv_chains

    pdb_code = row["pdb_code"].strip().upper()
    return PDB_CHAIN_OVERRIDES.get(
        pdb_code,
        DEFAULT_ALLOWED_CHAINS,
    )


def get_starting_chain(row: dict) -> str:
    """
    Priority:
      1. optional CSV column 'starting_chain'
      2. PDB-specific override (4CHA -> C)
      3. default chain A
    """
    csv_chain = (row.get("starting_chain", "") or "").strip()
    if csv_chain:
        return csv_chain

    pdb_code = row["pdb_code"].strip().upper()
    return PDB_STARTING_CHAIN_OVERRIDES.get(
        pdb_code,
        DEFAULT_STARTING_CHAIN,
    )


# UTILITIES

def normalize_text(value: str) -> str:
    """
    Normalize strings for robust filename matching.
    Removes punctuation/spaces/underscores/hyphens and ignores case.
    """
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


def safe_folder_name(value: str) -> str:
    """
    Make a readable Linux-safe folder name.
    """
    value = value.strip().replace("/", "_").replace("\\", "_")
    value = re.sub(r"[*?\"<>|:]", "", value)
    value = re.sub(r"\s+", "_", value)
    value = re.sub(r"_+", "_", value)
    return value.strip("._") or "unnamed_protein"


def pdb_code_aliases(pdb_code: str) -> list[str]:
    """Treat AF3 and AlphaFold3 as equivalent filename labels."""
    code = str(pdb_code).strip()
    aliases = [code]

    if normalize_text(code) in {"af3", "alphafold3"}:
        aliases.extend(["AF3", "AlphaFold3"])

    result = []
    seen = set()
    for item in aliases:
        key = normalize_text(item)
        if key not in seen:
            seen.add(key)
            result.append(item)
    return result


def strip_4rev_suffix(stem: str) -> str:
    return re.sub(r"(?i)_?4rev$", "", stem)


def resolve_exact_filename(value: str) -> Path:
    """
    Resolve an optional exact pdb_filename supplied in the CSV.
    Accepts either:
        1CEX_FsC_4rev.pdb
    or:
        1CEX_FsC_4rev
    """
    name = value.strip()
    if not name:
        raise ValueError("Empty pdb_filename.")

    candidate = PDB_DIR / name
    if candidate.exists():
        return candidate

    if not name.lower().endswith(".pdb"):
        candidate = PDB_DIR / f"{name}.pdb"
        if candidate.exists():
            return candidate

    raise FileNotFoundError(
        f"CSV pdb_filename points to a file that does not exist: {name}"
    )


def find_matching_pdb(
    protein_name: str,
    pdb_code: str,
    available_pdbs: list[Path],
    exact_filename: str = "",
) -> Path:
    """
    Priority:
    1. If pdb_filename is present in the CSV -> use it exactly.
    2. Otherwise match either PDBCODE_PROTEIN or PROTEIN_PDBCODE.

    Matching ignores punctuation and case.
    """
    if exact_filename.strip():
        return resolve_exact_filename(exact_filename)

    expected = set()
    for code in pdb_code_aliases(pdb_code):
        expected.add(normalize_text(f"{code}{protein_name}"))
        expected.add(normalize_text(f"{protein_name}{code}"))

    matches = []
    for path in available_pdbs:
        base = strip_4rev_suffix(path.stem)
        if normalize_text(base) in expected:
            matches.append(path)

    if len(matches) == 1:
        return matches[0]

    if len(matches) == 0:
        raise FileNotFoundError(
            f"No *_4rev.pdb matched nome_proteina='{protein_name}' "
            f"and pdb_code='{pdb_code}'. "
            f"If the real filename uses a different short name, add a "
            f"'pdb_filename' column to the CSV and write the exact filename."
        )

    raise RuntimeError(
        "Ambiguous filename match. More than one PDB matches this row: "
        + ", ".join(path.name for path in matches)
        + ". Add the exact 'pdb_filename' in the CSV."
    )


def canonicalize_pdb(src: Path, dst: Path, allowed_chains: tuple[str, ...]):
    """
    Create a CAVER-ready PDB copy.

    AMBER residue names are changed only in fixed PDB residue-name columns.
    Original input is NEVER modified.

    Returns:
        converted_residues
        remaining_noncanonical_ATOM_residues
    """
    dst.parent.mkdir(parents=True, exist_ok=True)

    converted = set()
    out_lines = []

    with src.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")) and len(line) >= 27:
                old_name = line[17:20].strip().upper()

                if old_name in AMBER_TO_STANDARD:
                    new_name = AMBER_TO_STANDARD[old_name]

                    chain = line[21].strip() or "-"
                    resseq = line[22:26].strip()
                    icode = line[26].strip()
                    resid = f"{resseq}{icode}"

                    converted.add((chain, resid, old_name, new_name))

                    # PDB residue name is columns 18-20 -> Python [17:20].
                    line = line[:17] + f"{new_name:>3}" + line[20:]

            out_lines.append(line)

    with dst.open("w", encoding="utf-8") as handle:
        handle.writelines(out_lines)

    remaining_noncanonical = set()
    atom_chains = set()

    # Check ATOM records only.
    # HETATM entries may legitimately be waters, ions, ligands, etc.
    with dst.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("ATOM  ") and len(line) >= 27:
                resname = line[17:20].strip().upper()
                atom_chains.add(line[21])

                if resname not in CANONICAL_AA:
                    chain = line[21].strip() or "-"
                    resseq = line[22:26].strip()
                    icode = line[26].strip()
                    remaining_noncanonical.add(
                        (resname, chain, f"{resseq}{icode}")
                    )

    expected_chain_set = set(allowed_chains)

    if atom_chains != expected_chain_set:
        shown = ", ".join(
            repr(c if c.strip() else "<blank>")
            for c in sorted(atom_chains)
        )
        expected = "+".join(allowed_chains)

        raise RuntimeError(
            f"Prepared CAVER input has the wrong chain set. "
            f"Expected chains: {expected}. "
            f"Found chain column values: {shown}. "
            f"For the main dataset use chain A; for 4CHA use exactly A+B+C."
        )

    return converted, remaining_noncanonical


def find_terminal_catalytic_hydrogen(
    pdb_path: Path,
    residue_number: int,
    expected_residue_name: str,
    expected_chain: str,
) -> tuple[int, str, str, str, float, float, float]:
    """
    Find ONLY the chemically terminal hydrogen of the catalytic side chain:

        SER -> HG
        THR -> HG1
        CYS -> HG

    The atom coordinates are then used as:
        starting_point_coordinates X Y Z

    Returns:
        atom_serial
        chain
        atom_name
        residue_identifier
        x
        y
        z
    """
    residue_name = expected_residue_name.strip().upper()

    if residue_name not in TERMINAL_CATALYTIC_HYDROGEN:
        raise RuntimeError(
            f"No terminal catalytic hydrogen rule is defined for "
            f"{residue_name}. Supported catalytic residues are "
            f"SER, THR and CYS."
        )

    target_atom = TERMINAL_CATALYTIC_HYDROGEN[residue_name]

    matches = []
    residue_ids = set()
    atoms_seen = set()
    same_number_other_names = set()

    with pdb_path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line.startswith("ATOM  ") or len(line) < 54:
                continue

            try:
                current_resnum = int(line[22:26].strip())
            except ValueError:
                continue

            if current_resnum != int(residue_number):
                continue

            current_resname = line[17:20].strip().upper()
            chain = line[21].strip()
            icode = line[26].strip()

            # Search the catalytic residue only in the requested chain.
            if chain != expected_chain:
                continue

            residue_id = (
                f"{chain}:{current_resnum}{icode}"
                if chain
                else f"{current_resnum}{icode}"
            )

            if current_resname != residue_name:
                same_number_other_names.add(
                    (residue_id, current_resname)
                )
                continue

            residue_ids.add(residue_id)
            atom_name = line[12:16].strip().upper()
            atoms_seen.add(atom_name)

            if atom_name != target_atom:
                continue

            try:
                serial = int(line[6:11].strip())
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
            except ValueError as exc:
                raise RuntimeError(
                    f"Cannot parse {residue_id}:{target_atom} in "
                    f"{pdb_path.name}: {exc}"
                )

            matches.append(
                (
                    serial,
                    chain,
                    atom_name,
                    residue_id,
                    x,
                    y,
                    z,
                )
            )

    if len(residue_ids) > 1:
        raise RuntimeError(
            f"Catalytic residue {residue_name}{residue_number} occurs in "
            f"multiple chains/residue IDs: {', '.join(sorted(residue_ids))}. "
            f"The script will not guess."
        )

    if len(matches) == 1:
        serial, chain, atom_name, residue_id, x, y, z = matches[0]

        if chain != expected_chain:
            raise RuntimeError(
                f"Invalid catalytic-atom chain for "
                f"{residue_id}:{atom_name}: found {chain!r}, "
                f"expected {expected_chain!r}."
            )

        return matches[0]

    if len(matches) == 0:
        detail = ""

        if residue_ids:
            detail = (
                f" Atoms found for {residue_name}{residue_number}: "
                + ", ".join(sorted(atoms_seen))
                + "."
            )
        elif same_number_other_names:
            detail = (
                " Residues with this number found instead: "
                + ", ".join(
                    f"{identifier}={name}"
                    for identifier, name
                    in sorted(same_number_other_names)
                )
                + "."
            )

        raise RuntimeError(
            f"Required terminal catalytic hydrogen {target_atom} of "
            f"{residue_name}{residue_number} was not found in "
            f"{pdb_path.name}.{detail} "
            f"No OG/OG1/SG fallback is allowed."
        )

    raise RuntimeError(
        f"Terminal catalytic hydrogen {target_atom} of "
        f"{residue_name}{residue_number} is ambiguous: "
        + ", ".join(
            f"{residue_id}:{atom_name} serial={serial}"
            for serial, chain, atom_name, residue_id, x, y, z
            in matches
        )
    )


def write_config(
    config_path: Path,
    x: float,
    y: float,
    z: float,
):
    """
    Write one traceable CAVER config per enzyme.

    The starting point is the terminal catalytic hydrogen:
        SER -> HG
        THR -> HG1
        CYS -> HG
        -> XYZ coordinates
        -> starting_point_coordinates X Y Z
    """
    config_path.parent.mkdir(parents=True, exist_ok=True)

    config = f"""# Automatically generated by run_caver_batch_last_atom_xyz.py
# Starting-point strategy: terminal catalytic hydrogen -> XYZ coordinates.
# SER->HG, THR->HG1, CYS->HG. No heavy-atom fallback.

# CALCULATION SETUP
load_tunnels no
load_cluster_tree no
stop_after never

# INPUT DATA
time_sparsity 1
first_frame 1
last_frame 100000

# TUNNEL CALCULATION
starting_point_coordinates {x:.15f} {y:.15f} {z:.15f}
probe_radius {PROBE_RADIUS}
shell_radius {SHELL_RADIUS}
shell_depth {SHELL_DEPTH}

# TUNNEL CLUSTERING
clustering average_link
weighting_coefficient 1
clustering_threshold 3.5
exclude_start_zone 2
exclude_end_zone 0
min_middle_zone 5
save_zones no

# GENERATION OF OUTPUTS
one_tunnel_in_snapshot cheapest
max_output_clusters 999
save_dynamics_visualization no
generate_summary yes
generate_tunnel_characteristics yes
generate_tunnel_profiles yes
generate_histograms no
compute_tunnel_residues yes
residue_contact_distance 3.0
compute_bottleneck_residues yes
bottleneck_contact_distance 3.0

# STARTING POINT OPTIMIZATION
max_distance {MAX_DISTANCE}
desired_radius {DESIRED_RADIUS}

# ADVANCED TUNNEL CALCULATION
number_of_approximating_balls 12
add_central_sphere yes
max_number_of_tunnels 10000
max_limiting_radius 100
cost_function_exponent 2
automatic_shell_radius no
automatic_shell_radius_bottleneck_multiplier 2
starting_point_protection_radius 4

# REDUNDANT TUNNELS REMOVAL
frame_clustering yes
frame_weighting_coefficient {FRAME_WEIGHTING_COEFFICIENT}
frame_clustering_threshold {FRAME_CLUSTERING_THRESHOLD}
frame_exclude_start_zone 0
frame_exclude_end_zone 0
frame_min_middle_zone 5

# AVERAGING OF TUNNEL ENDS
average_surface_frame yes
average_surface_global yes
average_surface_smoothness_angle 10
average_surface_point_min_angle 5
average_surface_tunnel_sampling_step 0.5

# APPROXIMATE CLUSTERING
do_approximate_clustering no
cluster_by_hierarchical_clustering 20000
max_training_clusters 15
generate_unclassified_cluster no

# OUTPUTS
profile_tunnel_sampling_step 0.5
visualization_tunnel_sampling_step 1
visualize_tunnels_per_cluster 5000
visualization_subsampling random
compute_errors no
save_error_profiles no
generate_trajectory no

# OTHERS
swap yes
seed 1
include_residue_names 20_AA
"""
    config_path.write_text(config, encoding="utf-8")


def summary_has_tunnels(output_dir: Path) -> bool:
    """
    Return True only if summary.txt contains at least one tunnel-cluster row.
    """
    summary_path = output_dir / "summary.txt"

    if not summary_path.exists():
        return False

    with summary_path.open(
        "r",
        encoding="utf-8",
        errors="replace",
    ) as handle:
        for line in handle:
            stripped = line.strip()

            # Tunnel-cluster rows begin with a numeric ID.
            if re.match(r"^\d+\s+", stripped):
                return True

    return False


def tunnels_stored_from_log(log_text: str) -> int:
    """
    Read CAVER messages such as:
        0 tunnels stored.
        18 tunnels stored.
    """
    counts = [
        int(value)
        for value in re.findall(
            r"(?im)^\s*(\d+)\s+tunnels?\s+stored\.",
            log_text,
        )
    ]

    return sum(counts)


def run_caver(
    input_dir: Path,
    config_path: Path,
    output_dir: Path,
    log_path: Path,
):

    if output_dir.exists():
        if OVERWRITE:
            shutil.rmtree(output_dir)
        else:
            raise FileExistsError(
                f"Output folder already exists: {output_dir}. "
                f"Set OVERWRITE=True only if you intentionally want to "
                f"delete and recompute it."
            )

    output_dir.parent.mkdir(parents=True, exist_ok=True)
    log_path.parent.mkdir(parents=True, exist_ok=True)

    command = [
        "java",
        f"-Xmx{JAVA_XMX}",
        "-cp", str(CAVER_HOME / "lib"),
        "-jar", str(CAVER_HOME / "caver.jar"),
        "-home", str(CAVER_HOME),
        "-pdb", str(input_dir),
        "-conf", str(config_path),
        "-out", str(output_dir),
    ]

    completed = subprocess.run(
        command,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )

    log_text = completed.stdout or ""
    log_path.write_text(log_text, encoding="utf-8")

    severe = "SEVERE:" in log_text
    starting_point_out = "Starting point is out" in log_text

    stored_tunnels = tunnels_stored_from_log(log_text)
    summary_tunnels = summary_has_tunnels(output_dir)
    tunnel_output = stored_tunnels > 0 or summary_tunnels

    return {
        "returncode": completed.returncode,
        "severe": severe,
        "starting_point_out": starting_point_out,
        "tunnel_output": tunnel_output,
        "stored_tunnels": stored_tunnels,
        "summary_has_tunnels": summary_tunnels,
        "command": " ".join(command),
    }


def load_csv_rows(csv_file: Path) -> list[dict]:
    required = {
        "nome_proteina",
        "pdb_code",
        "starting_residue_name",
        "starting_residue_number",
    }

    with csv_file.open(
        "r",
        newline="",
        encoding="utf-8-sig",
    ) as handle:
        reader = csv.DictReader(handle)

        if reader.fieldnames is None:
            raise RuntimeError("The CSV has no header.")

        missing = required - set(reader.fieldnames)
        if missing:
            raise RuntimeError(
                "The CSV is missing required columns: "
                + ", ".join(sorted(missing))
            )

        return list(reader)


def make_run_labels(rows: list[dict]) -> dict[int, str]:
    """
    Normally output folder = protein name.

    If a protein name occurs more than once, append pdb_code so that
    duplicated proteins do not overwrite each other.
    """
    counts = Counter(row["nome_proteina"].strip() for row in rows)
    labels = {}

    for index, row in enumerate(rows):
        protein = row["nome_proteina"].strip()
        pdb_code = row["pdb_code"].strip()

        if counts[protein] > 1:
            label = f"{protein}_{pdb_code}"
        else:
            label = protein

        labels[index] = safe_folder_name(label)

    return labels


def write_summary_csv(path: Path, rows: list[dict]):
    if not rows:
        # Still create a readable empty file.
        path.write_text("", encoding="utf-8")
        return

    with path.open(
        "w",
        newline="",
        encoding="utf-8-sig",
    ) as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(rows[0].keys()),
        )
        writer.writeheader()
        writer.writerows(rows)


# MAIN


def main():
    # -------------------------------------------------------------------------
    # Validate installation/input paths before processing anything.
    # -------------------------------------------------------------------------
    if not CSV_FILE.exists():
        raise FileNotFoundError(f"CSV not found: {CSV_FILE}")

    if not PDB_DIR.exists():
        raise FileNotFoundError(f"PDB directory not found: {PDB_DIR}")

    if not (CAVER_HOME / "caver.jar").exists():
        raise FileNotFoundError(
            f"CAVER jar not found: {CAVER_HOME / 'caver.jar'}"
        )

    if not (CAVER_HOME / "lib").exists():
        raise FileNotFoundError(
            f"CAVER lib directory not found: {CAVER_HOME / 'lib'}"
        )

    CAVER_TERMINAL_ROOT.mkdir(parents=True, exist_ok=True)

    # Technical pipeline files are separated from the actual per-enzyme
    # CAVER output directories.
    input_root = CAVER_TERMINAL_ROOT / "_pipeline_inputs"
    config_root = CAVER_TERMINAL_ROOT / "_pipeline_configs"
    log_root = CAVER_TERMINAL_ROOT / "_pipeline_logs"

    for path in (input_root, config_root, log_root):
        path.mkdir(parents=True, exist_ok=True)

    rows = load_csv_rows(CSV_FILE)
    run_labels = make_run_labels(rows)

    available_pdbs = sorted(PDB_DIR.glob("*_4rev.pdb"))
    if not available_pdbs:
        raise RuntimeError(
            f"No *_4rev.pdb files were found in {PDB_DIR}"
        )

    audit_rows = []

    # -------------------------------------------------------------------------
    # Process one enzyme at a time.
    # -------------------------------------------------------------------------
    for index, row in enumerate(rows):
        protein = row["nome_proteina"].strip()
        pdb_code = row["pdb_code"].strip()
        run_label = run_labels[index]

        allowed_chains = get_allowed_chains(row)
        starting_chain_expected = get_starting_chain(row)

        if ONLY_PROTEINS is not None and protein not in ONLY_PROTEINS:
            continue

        print()
        print("=" * 80)
        print(f"PROTEIN: {protein}")
        print(f"PDB/SOURCE: {pdb_code}")
        print(f"CAVER OUTPUT FOLDER: {run_label}")
        print(f"SELECTED CHAINS: {'+'.join(allowed_chains)}")
        print(f"CATALYTIC STARTING CHAIN: {starting_chain_expected}")
        print("=" * 80)

        audit = {
            "nome_proteina": protein,
            "pdb_code": pdb_code,
            "run_folder": run_label,
            "selected_chains": "+".join(allowed_chains),
            "starting_chain_expected": starting_chain_expected,
            "input_pdb": "",
            "prepared_pdb": "",
            "starting_residue_expected": "",
            "starting_atom_name": "",
            "starting_atom_serial": "",
            "starting_atom_chain": "",
            "starting_point_x": "",
            "starting_point_y": "",
            "starting_point_z": "",
            "starting_point_caver": "",
            "amber_residues_converted": "",
            "remaining_noncanonical_atom_residues": "",
            "caver_returncode": "",
            "caver_severe": "",
            "starting_point_outside": "",
            "tunnel_output_found": "",
            "tunnels_stored": "",
            "summary_has_tunnels": "",
            "status": "",
            "reason": "",
        }

        try:
            # -----------------------------------------------------------------
            # Read catalytic starting residue from CSV.
            # -----------------------------------------------------------------
            expected_resname = (
                row["starting_residue_name"].strip().upper()
            )

            raw_number = row["starting_residue_number"].strip()
            if not raw_number:
                raise ValueError(
                    "starting_residue_number is empty in the CSV."
                )

            residue_number = int(float(raw_number))

            audit["starting_residue_expected"] = (
                f"{expected_resname}{residue_number}"
            )

            # -----------------------------------------------------------------
            # Locate corresponding PDB.
            # -----------------------------------------------------------------
            exact_filename = row.get("pdb_filename", "") or ""

            input_pdb = find_matching_pdb(
                protein_name=protein,
                pdb_code=pdb_code,
                available_pdbs=available_pdbs,
                exact_filename=exact_filename,
            )

            audit["input_pdb"] = str(input_pdb)
            print(f"Input PDB: {input_pdb.name}")

            # -----------------------------------------------------------------
            # Prepare CAVER input without touching original PDB.
            # -----------------------------------------------------------------
            enzyme_input_dir = input_root / run_label

            if enzyme_input_dir.exists() and OVERWRITE:
                shutil.rmtree(enzyme_input_dir)

            enzyme_input_dir.mkdir(parents=True, exist_ok=True)

            # Keep the exact original filename for reproducibility.
            prepared_pdb = enzyme_input_dir / input_pdb.name

            converted, remaining_noncanonical = canonicalize_pdb(
                input_pdb,
                prepared_pdb,
                allowed_chains,
            )

            audit["prepared_pdb"] = str(prepared_pdb)

            converted_text = "; ".join(
                f"{chain}:{resid} {old}->{new}"
                for chain, resid, old, new in sorted(converted)
            )
            audit["amber_residues_converted"] = converted_text

            if converted:
                print(
                    f"AMBER residue conversion: "
                    f"{len(converted)} residue(s)"
                )
                for chain, resid, old, new in sorted(converted):
                    print(f"  {chain}:{resid}  {old} -> {new}")
            else:
                print("AMBER residue conversion: none required")

            remaining_text = "; ".join(
                f"{name} {chain}:{resid}"
                for name, chain, resid
                in sorted(remaining_noncanonical)
            )
            audit[
                "remaining_noncanonical_atom_residues"
            ] = remaining_text

            if remaining_noncanonical:
                raise RuntimeError(
                    "Non-canonical ATOM residues remain after the known "
                    "AMBER conversions: "
                    + remaining_text
                    + ". The enzyme was NOT sent to CAVER."
                )

            # -----------------------------------------------------------------
            # Use ONLY the terminal hydrogen of the catalytic side chain:
            # SER -> HG, THR -> HG1, CYS -> HG.
            # No heavy-atom fallback.
            # -----------------------------------------------------------------
            (
                starting_serial,
                starting_chain,
                starting_atom_name,
                starting_residue_id,
                starting_x,
                starting_y,
                starting_z,
            ) = find_terminal_catalytic_hydrogen(
                prepared_pdb,
                residue_number,
                expected_resname,
                starting_chain_expected,
            )

            audit["starting_atom_name"] = starting_atom_name
            audit["starting_atom_serial"] = starting_serial
            audit["starting_atom_chain"] = starting_chain or "-"

            chain_for_print = starting_chain or "-"

            print(
                f"Starting atom (terminal catalytic hydrogen): "
                f"{chain_for_print}:{residue_number}:{starting_atom_name} "
                f"(serial {starting_serial})"
            )

            # -------------------------------------------------------------
            # Manual starting-point override ONLY for selected proteins.
            # The catalytic atom is still identified and reported above,
            # but CAVER uses the manually chosen XYZ coordinates.
            # -------------------------------------------------------------
            if protein in MANUAL_STARTING_POINTS:
                starting_x, starting_y, starting_z = (
                    MANUAL_STARTING_POINTS[protein]
                )

                print(
                    f"MANUAL starting-point override for {protein}: "
                    f"{starting_x:.6f} "
                    f"{starting_y:.6f} "
                    f"{starting_z:.6f}"
                )
            else:
                print(
                    f"Starting XYZ from terminal catalytic hydrogen: "
                    f"{starting_x:.3f} "
                    f"{starting_y:.3f} "
                    f"{starting_z:.3f}"
                )

            audit["starting_point_x"] = f"{starting_x:.15f}"
            audit["starting_point_y"] = f"{starting_y:.15f}"
            audit["starting_point_z"] = f"{starting_z:.15f}"
            audit["starting_point_caver"] = (
                f"starting_point_coordinates "
                f"{starting_x:.15f} "
                f"{starting_y:.15f} "
                f"{starting_z:.15f}"
            )

            # -----------------------------------------------------------------
            # Write one config file per enzyme.
            # -----------------------------------------------------------------
            config_path = (
                config_root / f"{run_label}_config.txt"
            )
            write_config(
                config_path,
                starting_x,
                starting_y,
                starting_z,
            )

            print(f"Config: {config_path}")

           
            # The actual CAVER output is named by protein, NOT by run number.
         
            output_dir = CAVER_TERMINAL_ROOT / run_label
            log_path = log_root / f"{run_label}.log"

            print(f"Running CAVER -> {output_dir}")

            result = run_caver(
                input_dir=enzyme_input_dir,
                config_path=config_path,
                output_dir=output_dir,
                log_path=log_path,
            )

            audit["caver_returncode"] = result["returncode"]
            audit["caver_severe"] = result["severe"]
            audit[
                "starting_point_outside"
            ] = result["starting_point_out"]
            audit[
                "tunnel_output_found"
            ] = result["tunnel_output"]
            audit["tunnels_stored"] = result["stored_tunnels"]
            audit["summary_has_tunnels"] = result[
                "summary_has_tunnels"
            ]

            # -----------------------------------------------------------------
            # Conservative success criterion.
            # -----------------------------------------------------------------
            if (
                result["returncode"] == 0
                and not result["severe"]
                and not result["starting_point_out"]
                and result["tunnel_output"]
            ):
                audit["status"] = "OK"
                audit["reason"] = ""
                print("STATUS: OK")

            else:
                audit["status"] = "FAILED"

                reasons = []

                if result["returncode"] != 0:
                    reasons.append(
                        f"Java/CAVER return code "
                        f"{result['returncode']}"
                    )

                if result["severe"]:
                    reasons.append(
                        "CAVER log contains SEVERE"
                    )

                if result["starting_point_out"]:
                    reasons.append(
                        "CAVER reports that the starting point is outside"
                    )

                if not result["tunnel_output"]:
                    reasons.append(
                        "CAVER stored 0 tunnels / summary.txt contains "
                        "no tunnel cluster"
                    )

                audit["reason"] = (
                    "; ".join(reasons)
                    + f". See log: {log_path}"
                )

                print("STATUS: FAILED")
                print(f"REASON: {audit['reason']}")

        except FileExistsError as exc:
            audit["status"] = "SKIPPED"
            audit["reason"] = str(exc)

            print("STATUS: SKIPPED")
            print(f"REASON: {exc}")

        except Exception as exc:
            audit["status"] = "FAILED"
            audit["reason"] = str(exc)

            print("STATUS: FAILED")
            print(f"REASON: {exc}")

        audit_rows.append(audit)


    # Final audit files.

    audit_path = CAVER_TERMINAL_ROOT / "caver_batch_audit.csv"

    successful_rows = [
        row for row in audit_rows
        if row["status"] == "OK"
    ]

    failed_rows = [
        row for row in audit_rows
        if row["status"] == "FAILED"
    ]

    skipped_rows = [
        row for row in audit_rows
        if row["status"] == "SKIPPED"
    ]

    failed_path = (
        CAVER_TERMINAL_ROOT / "caver_failed_proteins.csv"
    )
    successful_path = (
        CAVER_TERMINAL_ROOT / "caver_successful_proteins.csv"
    )

    write_summary_csv(audit_path, audit_rows)
    write_summary_csv(failed_path, failed_rows)
    write_summary_csv(successful_path, successful_rows)


    print()
    print("#" * 80)
    print("CAVER BATCH SUMMARY")
    print("#" * 80)
    print(f"Processed:  {len(audit_rows)}")
    print(f"Successful: {len(successful_rows)}")
    print(f"Failed:     {len(failed_rows)}")
    print(f"Skipped:    {len(skipped_rows)}")

    if failed_rows:
        print()
        print("PROTEINS THAT COULD NOT BE PROCESSED:")
        print("-" * 80)

        for row in failed_rows:
            print(
                f"- {row['nome_proteina']} "
                f"({row['pdb_code']}): {row['reason']}"
            )
    else:
        print()
        print("No failed proteins.")

    if skipped_rows:
        print()
        print("SKIPPED PROTEINS:")
        print("-" * 80)

        for row in skipped_rows:
            print(
                f"- {row['nome_proteina']} "
                f"({row['pdb_code']}): {row['reason']}"
            )

    print()
    print(f"Main output directory: {CAVER_TERMINAL_ROOT}")
    print(f"Full audit:            {audit_path}")
    print(f"Failed proteins:       {failed_path}")
    print(f"Successful proteins:   {successful_path}")
    print("#" * 80)


if __name__ == "__main__":
    main()
