import argparse
import os
from pathlib import Path
from shutil import copyfile, which
import subprocess
from typing import Union
import glob

from surfmap.bin.multival_csv_to_pdb import run as run_multival_csv2pdb
from surfmap.lib.logs import get_logger


logger = get_logger(name=__name__)

def _estimate_pqr_grid_lengths(pqr_file: Union[str, Path], extra_padding: float=20.0):
    """Estimate APBS coarse/fine grid lengths from PQR coordinates.

    Returns:
        tuple[float, float, float]: lengths for x, y, z dimensions.
    """
    xs, ys, zs = [], [], []

    with open(pqr_file) as _f:
        for line in _f:
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            cols = line.split()
            if len(cols) < 6:
                continue

            # PQR (whitespace) expected x,y,z around columns 6,7,8.
            # Keep a conservative parser with fallback to PDB-like fixed columns.
            try:
                x, y, z = float(cols[5]), float(cols[6]), float(cols[7])
            except Exception:
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                except Exception:
                    continue

            xs.append(x)
            ys.append(y)
            zs.append(z)

    if not xs:
        logger.error(f"Could not parse coordinates from generated PQR file: {pqr_file}")
        exit(1)

    lx = (max(xs) - min(xs)) + extra_padding
    ly = (max(ys) - min(ys)) + extra_padding
    lz = (max(zs) - min(zs)) + extra_padding
    return max(lx, 20.0), max(ly, 20.0), max(lz, 20.0)


def _write_fallback_apbs_input(apbs_input_file: Union[str, Path], pqr_file: Union[str, Path], pqr_stem: str):
    """Write a minimal APBS input file when inputgen.py fails."""
    lx, ly, lz = _estimate_pqr_grid_lengths(pqr_file=pqr_file)

    content = f"""read
    mol pqr {pqr_file}
end
elec
    mg-auto
    dime 161 161 161
    cglen {lx:.3f} {ly:.3f} {lz:.3f}
    fglen {max(lx*0.6, 15.0):.3f} {max(ly*0.6, 15.0):.3f} {max(lz*0.6, 15.0):.3f}
    cgcent mol 1
    fgcent mol 1
    mol 1
    lpbe
    bcfl sdh
    pdie 2.0
    sdie 78.54
    chgm spl2
    srfm smol
    srad 1.4
    swin 0.3
    temp 298.15
    calcenergy no
    calcforce no
    write pot dx {pqr_stem}
end
quit
"""

    with open(apbs_input_file, "w") as _out:
        _out.write(content)


def _resolve_apbs_tools():
    """Resolve APBS executables/scripts paths from env var APBS or PATH.

    Returns:
        tuple[str, str, str]: paths for apbs, inputgen.py, multivalue
    """
    apbs_env = os.getenv("APBS")

    if apbs_env:
        path_apbs = Path(apbs_env)
        apbs = str(path_apbs / "bin" / "apbs")
        inputgen = str(path_apbs / "share" / "apbs" / "tools" / "manip" / "inputgen.py")
        multivalue = str(path_apbs / "share" / "apbs" / "tools" / "bin" / "multivalue")

        if Path(apbs).exists() and Path(inputgen).exists() and Path(multivalue).exists():
            return apbs, inputgen, multivalue

        logger.error(
            "Environment variable APBS is set, but required files were not found under it. "
            "Expected: bin/apbs, share/apbs/tools/manip/inputgen.py, share/apbs/tools/bin/multivalue"
        )
        exit(1)

    apbs_bin = which("apbs")
    if not apbs_bin:
        logger.error("Could not find 'apbs' executable. Please install APBS or set APBS env variable to APBS root path.")
        exit(1)

    apbs_path = Path(apbs_bin).resolve()
    apbs_root = apbs_path.parent.parent

    # First, try APBS-relative layout (works for many package/install variants)
    inputgen_candidates = [
        apbs_root / "share" / "apbs" / "tools" / "manip" / "inputgen.py",
        apbs_root / "share" / "apbs" / "tools" / "inputgen.py",
    ]
    multivalue_candidates = [
        apbs_root / "share" / "apbs" / "tools" / "bin" / "multivalue",
        apbs_root / "share" / "apbs" / "tools" / "multivalue",
    ]

    # Then try executables/scripts that may already be directly exposed in PATH
    inputgen_in_path = which("inputgen.py")
    multivalue_in_path = which("multivalue")
    if inputgen_in_path:
        inputgen_candidates.append(Path(inputgen_in_path))
    if multivalue_in_path:
        multivalue_candidates.append(Path(multivalue_in_path))

    # Lastly, search common system prefixes (apt/conda/module setups can vary)
    for pattern in [
        "/usr/**/apbs/**/inputgen.py",
        "/opt/**/apbs/**/inputgen.py",
        str(apbs_root / "**" / "inputgen.py"),
    ]:
        inputgen_candidates.extend(Path(p) for p in glob.glob(pattern, recursive=True))

    for pattern in [
        "/usr/**/apbs/**/multivalue",
        "/opt/**/apbs/**/multivalue",
        str(apbs_root / "**" / "multivalue"),
    ]:
        multivalue_candidates.extend(Path(p) for p in glob.glob(pattern, recursive=True))

    inputgen = next((str(p) for p in inputgen_candidates if p.exists()), None)
    multivalue = next((str(p) for p in multivalue_candidates if p.exists()), None)

    # if not inputgen or not multivalue:
    #     logger.error(
    #         "Could not resolve APBS helper tools (inputgen.py/multivalue). "
    #         f"Detected apbs binary: {apbs_path}. "
    #         "Please set APBS env variable to APBS root path containing share/apbs/tools."
    #     )
    #     exit(1)

    return str(apbs_path), inputgen, multivalue

def get_args():
    parser = argparse.ArgumentParser(
        description="Script to compute electrostatics. Returns PDB files with coordinates of the shell surface particles and electrostatic.")

    parser.add_argument(
        "-pdb",
        required=True,
        type=str,
        help="Path to the pdbfile of the protein complex"
    )

    parser.add_argument(
        "-csv",
        required=True,
        type=str,
        help="CSV file with coordinates of the shell particles."
    )

    parser.add_argument(
        "-pqr",
        required=False,
        type=str,
        default="",
        help="PQR file used to build input file of APBS"
    )

    parser.add_argument(
        "-out",
        type=str,
        default=Path("."),
        help="Path of the output"
    )

    parser.add_argument(
        "-subdir",
        type=str,
        default=Path("tmp-elec"),
        help="Path name of the electrostatics directory content"
    )

    return parser.parse_args()


def edit_inputgen(inputgen_file: Union[str, Path], pqrfile: Union[str, Path]):
    with open(inputgen_file, 'r') as _file:
        data = _file.read()

    # Replace some data
    data = data.replace(Path(pqrfile).name, pqrfile)
    data = data.replace("pot dx pot", f"pot dx {pqrfile}")

    # Write the file out
    with open(inputgen_file, 'w') as _file:
        _file.write(data)

def _resolve_potential_grid_file(outdir: Union[str, Path], outfile_basename: str, outfile_pqr: str):
    outdir = Path(outdir)

    candidates = [
        Path(f"{outfile_pqr}.dx"),
        Path(f"{outfile_basename}.dx"),
    ]

    for candidate in candidates:
        if candidate.exists():
            return str(candidate)

    dx_files = sorted(outdir.glob("*.dx"), key=lambda p: p.stat().st_mtime, reverse=True)
    if dx_files:
        logger.warning(
            f"Expected APBS grid file not found under default names. Using detected file: {dx_files[0]}"
        )
        return str(dx_files[0])

    return None



def run(pdb_filename: Union[str, Path], csv_filename: Union[str, Path], force_field: str = "CHARMM",
        pqr_filename: Union[str, Path] = "", out_dir: Union[str, Path] = "."):
    """Calls the function run_compute_electrostatics() to compute electrostatics potential of the PDB file.

    Args:
        pdb_filename (Union[str, Path]): Path to a PDB filename
        csv_filename (Union[str, Path]): Path to a CSV filename (output of run_compute_shell function)
        pqr_filename (Union[str, Path]): Optional path to a PQR filename
        out_dir (Union[str, Path]): Path to the output directory

    Returns:
        str: Path to a PDB file with coordinates of the shell particles and atomic electrostatic potential values in bfactor
    """
    PDB2PQR_FORCE_FIELDS = ["AMBER", "CHARMM", "PARSE", "TYL06", "PEOEPB", "SWANSON"]

    # check if force field is allowed
    if force_field.upper() not in PDB2PQR_FORCE_FIELDS:
        logger.error(f"Error, the force field {force_field} is not accepted.")
        logger.error(f"Accepted force fields are: {' '.join(PDB2PQR_FORCE_FIELDS)}.")
        logger.error(f"Exiting now.\n")
        exit(1)

    # define access to useful APBS executables
    # PATH_APBS = Path(os.getenv("APBS"))
    # apbs = f"{str(PATH_APBS / 'bin' / 'apbs')}"
    # inputgen = f"{str(PATH_APBS)}/share/apbs/tools/manip/inputgen.py"
    # multivalue = f"{str(PATH_APBS)}/share/apbs/tools/bin/multivalue"
    apbs, inputgen, multivalue = _resolve_apbs_tools()

    # Create directory containing electrostatics calculation related files
    outdir = Path(out_dir)
    outdir.mkdir(exist_ok=True, parents=True)
    outfile_basename = str(outdir / Path(pdb_filename).stem)

    # compute pqr file if not given as input
    outfile_pqr = f"{outfile_basename}.pqr"
    if not pqr_filename:
        cmd_pdb2pqr = ["pdb2pqr30", "--ff", force_field.upper(), "--whitespace", pdb_filename, outfile_pqr]
        logger.debug("Convert PDB to PQR format")
        status = subprocess.run(cmd_pdb2pqr, capture_output=True)
        if status.returncode != 0:
            logger.error(
                f"Error occured during the conversion of PDB to PQR, the process will stop. Status: {status.returncode}")
            exit(1)
    else:
        logger.debug(f"Copying user-given PQR file as {outfile_pqr}")
        copyfile(src=pqr_filename, dst=outfile_pqr)


    # generate APBS input file
    outfile_in = f"{outfile_basename}.in"  # automaticall generated in the same path as the pqr file
    cmd_inputgen = ["python3", inputgen, "--method=auto", outfile_pqr]
    logger.debug(f"Running inputgen command: {' '.join(cmd_inputgen)}")
    status = subprocess.run(args=cmd_inputgen, capture_output=True, text=True)
    if status.returncode != 0:
        logger.warning(
            f"inputgen command failed (status: {status.returncode}). Falling back to internal APBS input writer.")
        if status.stderr:
            logger.warning(f"inputgen stderr: {status.stderr.strip()}")
        if status.stdout:
            logger.warning(f"inputgen stdout: {status.stdout.strip()}")
        logger.debug(f"Writing fallback APBS input file: {outfile_in}")
        _write_fallback_apbs_input(apbs_input_file=outfile_in, pqr_file=outfile_pqr, pqr_stem=outfile_pqr)
    else:
        logger.debug(f"Editing generated inputgen file: {outfile_in}")
        edit_inputgen(inputgen_file=outfile_in, pqrfile=outfile_pqr)

    # compute APBS electrostatics calculation
    outfile_pot = f"{outfile_pqr}.dx"  # default output name of APBS
    cmd_apbs = [apbs, outfile_in]
    logger.debug(f"Running APBS command: {cmd_apbs}")
    status = subprocess.run(args=cmd_apbs, capture_output=True, text=True)
    if status.returncode != 0:
        logger.error(f"Error occured during the APBS command, the process will stop. Status: {status.returncode}")
        if status.stderr:
            logger.error(f"APBS stderr: {status.stderr.strip()}")
        if status.stdout:
            logger.error(f"APBS stdout: {status.stdout.strip()}")

        exit(1)

    # cross the MSMS generated CSV file with the potential gridfile created by APBS
    outfile_multivalue = f"{outfile_basename}.mult"
    if not Path(csv_filename).exists():
        logger.error(f"Could not find shell coordinates file required by multivalue: {csv_filename}")
        exit(1)
    resolved_pot = _resolve_potential_grid_file(outdir=outdir, outfile_basename=outfile_basename,
                                                outfile_pqr=outfile_pqr)
    if not resolved_pot:
        logger.error(
            "Could not find APBS potential grid file (*.dx) required by multivalue in output directory: "
            f"{outdir}"
        )
        exit(1)

    cmd_multivalue = [multivalue, csv_filename, resolved_pot, outfile_multivalue]
    if not os.access(multivalue, os.X_OK):
        cmd_multivalue = ["python3", multivalue, csv_filename, resolved_pot, outfile_multivalue]

    logger.debug(f"Running multivalue command: {cmd_multivalue}")
    status = subprocess.run(args=cmd_multivalue, capture_output=True, text=True)
    if status.returncode != 0:
        logger.error(f"Error occured during the multivalue command, the process will stop. Status: {status.returncode}")
        if status.stderr:
            logger.error(f"multivalue stderr: {status.stderr.strip()}")
        if status.stdout:
            logger.error(f"multivalue stdout: {status.stdout.strip()}")
        logger.error(
            "Hints: verify APBS/multivalue compatibility and that APBS generated a valid .dx file for your PDB. "
            "If needed, set APBS env variable to the exact APBS root used to run 'apbs'."
        )
        exit(1)

    # convert .csv file into a PDB file
    outfile_pdb = f"{outfile_basename}_shell.pdb"
    logger.debug(f"Converting CSV to PDB file")
    run_multival_csv2pdb(ifile=outfile_multivalue, ofile=outfile_pdb)

    # remove io.mc generated by APBS
    io_apbs = Path.cwd() / "io.mc"
    logger.debug(f"Removing APBS derived log file: {io_apbs}")
    io_apbs.unlink(missing_ok=True)

    return outfile_pdb


def main():
    args = get_args()

    return run(
        pdb_filename=args.pdb,
        csv_filename=args.csv,
        pqr_filename=args.pqr,
        out_dir=args.out,
        out_subdir=args.subdir
    )

