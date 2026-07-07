#!/usr/bin/env python
# -*- coding: utf-8 -*-
import argparse
from pathlib import Path
import sys
from typing import Union

from surfmap.tools import Structure
from surfmap.lib.logs import get_logger

logger = get_logger(name=__name__)

MASKED_ELECTROSTATICS_VALUE = -9999.0


def load_residue_mask(maskfile: Union[str, Path]):
    """Load residue identifiers to mask from a text file.

    Expected format: one residue per line with 3 columns separated by spaces/tabs:
    chainid resid resname
    """
    residues_to_mask = set()

    with open(maskfile) as maskf:
        for nline, line in enumerate(maskf, start=1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue

            cols = stripped.split()
            if len(cols) < 3:
                logger.error(
                    f"Error in electrostatics mask file on line {nline}: expected 3 columns (chainid resid resname)."
                )
                sys.exit(1)

            chain, resid, resname = cols[:3]
            residues_to_mask.add((chain, resid, resname.upper()))

    return residues_to_mask

def _replace_pdb_bfactor(line: str, value: float) -> str:
    """Return a PDB line with only occupancy/B-factor columns replaced."""
    line_body = line.rstrip("\n")
    padded = line_body.ljust(67)
    return f"{padded[:54]}{1.00:6.2f}{value:7.3f}{padded[67:]}\n"


def write_shell_pdb(dshell: dict, outfile: str, source_shell: Union[str, Path]=None):
    """Write shell particles as a PDB with mapped values in the B-factor column.

    If source_shell is provided, its ATOM/HETATM lines are copied and only the
    occupancy/B-factor columns are replaced. This keeps the original shell
    geometry, atom order, serials and formatting intact.
    """
    if source_shell:
        particles = iter(dshell["partlist"])
        written_particles = 0
        with open(source_shell) as in_shell, open(outfile, "w") as out_shell:
            for line in in_shell:
                if line.startswith(("ATOM", "HETATM")):
                    try:
                        particle = next(particles)
                    except StopIteration:
                        logger.error("Source shell contains more ATOM/HETATM records than mapped shell particles.")
                        sys.exit(1)
                    value = float(dshell[particle].get("charge", 0.0))
                    out_shell.write(_replace_pdb_bfactor(line, value))
                    written_particles += 1
                else:
                    out_shell.write(line)

            try:
                next(particles)
            except StopIteration:
                pass
            else:
                logger.error("Mapped shell contains more particles than source shell ATOM/HETATM records.")
                sys.exit(1)

        logger.debug(f"Wrote {written_particles} mapped shell particles while preserving source shell geometry")
        return

# def write_shell_pdb(dshell: dict, outfile: str):
#     with open(outfile, "w") as out_shell:
#         for i, particle in enumerate(dshell["partlist"], start=1):
#             # line = "{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}\n".format(
#             #     "ATOM",
#             #     i,
#             #     "O",
#             #     "",
#             #     "POS",
#             #     "A",
#             #     1,
#             #     "",
#             #     dshell[particle]["x"],
#             #     dshell[particle]["y"],
#             #     dshell[particle]["z"],
#             #     1.00,
#             #     dshell[particle]["charge"]
#             # )
#             value = float(dshell[particle].get("charge", 0.0))
#             line = (
#                 f"ATOM  {i:5d}  O   POS A   1    "
#                 f"{dshell[particle]['x']:8.3f}{dshell[particle]['y']:8.3f}{dshell[particle]['z']:8.3f}"
#                 f"{1.00:6.2f}{value:7.3f}          O\n"
#             )
#             out_shell.write(line)
#         out_shell.write("END\n")

def _copy_shell_values_by_order(target_shell: dict, value_shell: dict):
    """Copy per-particle values while preserving target shell geometry/order."""
    if len(target_shell["partlist"]) != len(value_shell["partlist"]):
        logger.error(
            "Cannot combine shell coordinates and shell values: particle counts differ "
            f"({len(target_shell['partlist'])} vs {len(value_shell['partlist'])})."
        )
        sys.exit(1)

    for target_particle, value_particle in zip(target_shell["partlist"], value_shell["partlist"]):
        target_shell[target_particle]["charge"] = value_shell[value_particle]["charge"]

def run_particles_mapping(shell: Union[str, Path], pdb: Union[str, Path], tomap: str, outdir: Union[str, Path]=".", res: Union[str, Path]=None, electrostatics_mask: Union[str, Path]=None, reference_shell: Union[str, Path]=None, write_mapped_shell_pdb: bool=False):
    """Assign atom/residue property value to its closest shell particle and compute the spherical coordinates of this particle.

    Two PDB structures are required:
        - one with coordinates of the shell particles
        - one with coordinates of the initial atoms of a protein

    Each shell particle must be assigned a property value. With the exception of the electrostatics property, each value assigned to a particle is retrieved from the property value
    of its closest atom/residue of the protein. For electrostatics property, values are already assigned to each particle.

    For hydrophobic properties (stickiness, kyte_doolittle and wimley_white), the values are computed on the fly for the closest atom/residue of a given particle.
    For properties like interface, circular_variance or circular_variance_atom, bfactor, the values are retrieved in the bfactor column. So, when circular variance property is asked,
    a new PDB is first generated.

    If a file with a list of residues to map is given, an output file with '_sph_coords.out' as suffix/extension is generated. This file could then be given to the
    compute_map function to map those residues.

    Args:
        shell (Union[str, Path]): path to the shell file produced with MSMS
        pdb (Union[str, Path]): Path to the input pdb file
        tomap (str): Property to map
        outdir (str): Path to the output directory. Defaults to '.'.
        res (Union[str, Path]): Optional path to a file with list of residues to map (optional)
        electrostatics_mask (Union[str, Path]): Optional path to a residue mask file used only for electrostatics.

    Returns:
        outfile_res_to_map (str | None): path to the spherical coords file of the list of residues to be mapped. None if no list of residues exists.
        partlist_outfile (str): Path to the file with spherical coords of every shell particles. The file has the following format: #phi #theta #scalevalue, #resnb, #resname, #chain

    """
    if not Path(shell).exists():
        print("Could not find the shell file. Exiting now\n")
        exit()

    # extracting information from receptor and particles
    is_bfactor = False if "electrostatics" in tomap else True
    is_charge = True if "electrostatics" in tomap else False

    if "circular_variance" in tomap:
        perres = False if "atom" in tomap else True
        pdb_cv_filename = str(Path(outdir) / f"{Path(pdb).stem}_{tomap}.pdb")
        logger.debug(f"Running circular variance computation of {pdb_cv_filename}")
        Structure.compute_CV(pdb, perres=perres, outfilename=pdb_cv_filename)

    # get pdb dict structure of input pdb and shell particles
    input_pdb = pdb_cv_filename if "circular_variance" in tomap else pdb

    logger.debug(f"Get dictionary of the PDB structure from {input_pdb}")
    dPDB = Structure.parsePDBMultiChains(input_pdb, bfactor=is_bfactor)

    logger.debug(f"Get dictionary of the shell structure from {shell}")
    dshell = Structure.parsePDBParticule(shell, infile=is_charge)

    if reference_shell and Path(reference_shell) != Path(shell):
        logger.debug(f"Preserving shell coordinates/order from reference shell {reference_shell}")
        value_shell = dshell
        dshell = Structure.parsePDBParticule(reference_shell, infile=False)
        _copy_shell_values_by_order(target_shell=dshell, value_shell=value_shell)

    # remove pdb of circular_variance if exists
    try:
        Path(pdb_cv_filename).unlink(missing_ok=True)
        logger.debug(f"{pdb_cv_filename} has been removed")
    except:
        pass

    # redefine tomap if necessary
    if tomap in ["interface", "circular_variance", "circular_variance_atom"]:
        property = "bfactor"
    else:
        property = tomap

    residues_to_mask = set()
    if electrostatics_mask and property == "electrostatics":
        logger.debug(f"Loading electrostatics residue mask from {electrostatics_mask}")
        residues_to_mask = load_residue_mask(electrostatics_mask)

    #  generate specific spherical coords file if some residues are asked to be mapped
    outfile_res_to_map = None
    if res:
        outfile_res_to_map = str(Path(outdir) / f"{Path(res).stem}_sph_coords.out")
        logger.debug(f"Parsing user-given residues to map and generating {outfile_res_to_map}")
        parseResidueList(res, dPDB, outfile=outfile_res_to_map)

    # compute center of mass of receptor
    logger.debug(f"Retrieving the center of mass of the pdb structure")
    CMR = Structure.centerMassResidueList(dPDB, all=True, reslist=False, computeForMulti=True, chain=" ")

    # get coords structure to compute distances
    coordlist, idres = Structure.get_coords_idres(dPDB=dPDB)

    partlist_outfile = Path(outdir) / f"{Path(pdb).stem}_{tomap}_partlist.out"
    with open(partlist_outfile, "w") as out:
        out.write(f"{'phi'}\t{'theta'}\t{'value'}\t{'resnb'}\t{'restype'}\t{'chain'}\n")

        # looping over all particules and computing corresponding Value
        logger.trace(
            f"Looping over all shell particles to assign them the property value of their closest atoms/residues")
        with open('closest_atoms_to_shell_pts.txt', 'w', encoding='utf-8') as out1:
            for i, particle in enumerate(dshell["partlist"], start=1):

                # get coords of the particule
                logger.trace(f"Reading coordinates of the particle {i} and converting it into spherical coordinates")
                coord_particle = (dshell[particle]["x"], dshell[particle]["y"], dshell[particle]["z"])
                _, phi, theta = Structure.coord2spherical(CMR, coord_particle)

                # get the closest atom of the shell particle
                logger.trace(f"Retrieving the closest atom of the particle {i}")
                _, closestAtom = Structure.getAtomCMRDist(coordlist=coordlist, idres=idres, CMR=coord_particle)
                chainid, resid, atomtype = closestAtom.split("_")
                # print(f"{closestAtom} for shell_id: {i}, {particle}")

                out1.write(f"{closestAtom} for shell_id: {i} \n")

                # assign the property value of atoms/residues to their closest shell particle or directly from the shell particle
                if property == "electrostatics":
                    scalevalue = dshell[particle]["charge"]
                    logger.trace(
                        f"Electrostatic value {scalevalue} read from the shell structure has been assigned to the particle {i}")
                elif property == "bfactor":
                    scalevalue = dPDB[chainid][resid][atomtype]["bfactor"]
                    logger.trace(
                        f"{tomap} value {scalevalue} read from the atom {chainid}-{resid}-{atomtype} of {input_pdb} has been assigned to the particle {i}")
                else:
                    scalevalue = Structure.returnPropensity(aa=dPDB[chainid][resid]["resname"], scale=tomap)
                    logger.trace(
                        f"{tomap} value {scalevalue} computed for the residue {chainid}-{resid} of {input_pdb} has been assigned to the particle {i}")

                resname = dPDB[chainid][resid]["resname"]
                resnum = dPDB[chainid][resid]["resnum"]

                is_masked = property == "electrostatics" and residues_to_mask and (chainid, str(resnum),
                                                                                   str(resname).upper()) in residues_to_mask
                if is_masked:
                    logger.trace(
                        f"Masking electrostatics value for residue {chainid}-{resnum}-{resname} "
                        f"by setting it to out-of-range sentinel {MASKED_ELECTROSTATICS_VALUE}"
                    )
                    scalevalue = MASKED_ELECTROSTATICS_VALUE

                # temporary - assign value to shell
                dshell[particle]["charge"] = scalevalue

                # out.write(
                #     "{:8f} {:8f} {:3f} {:8} {:8} {:8}\n".format(phi, theta, scalevalue, resnum, resname, chainid))
                out.write(f"{phi:8f} {theta:8f} {scalevalue:3f} {resnum:8} {resname:8} {chainid:8}\n")

    # temporary - write shell pdb with assigned values
    #shell_pdb = Path(outdir) / f"{Path(pdb).stem}_{tomap}_shell.pdb"
    #write_shell_pdb(dshell=dshell, outfile=shell_pdb)
    #print(dshell)
    if write_mapped_shell_pdb:
        shell_pdb = Path(outdir) / f"{Path(pdb).stem}_{tomap}_mapped_shell.pdb"
        logger.debug(f"Writing mapped shell PDB to {shell_pdb}")
        write_shell_pdb(dshell=dshell, outfile=shell_pdb, source_shell=reference_shell or shell)

    return outfile_res_to_map, partlist_outfile


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-pdb", required=True, help="input pdb file (path + file name)")
    parser.add_argument("-shell", required=True, help="shell file produced with MSMS")
    parser.add_argument("-d", required=True, default=".", help="output directory")
    parser.add_argument("-res", required=False, default=None, help="residue list to map (optional)")
    parser.add_argument("-tomap", required=False, type=str, default="electrostatics",
                        help="value to map (stickiness, kyte-doolittle, wimley-white, bfactor or elec. (Default: elec).")
    parser.add_argument("-surfth", required=False, type=float, default="5.0",
                        help="output file name (default = map+option in input directory)")
    parser.add_argument("-dist", required=False,
                        help="distance threshold that separates the spheres from the protein surface (default 5A)")
    parser.add_argument("-onlySurf", required=False, action="store_true",
                        help="to use only surface residues for the calculation of the value to map (Default: False).")
    args = parser.parse_args()

    shell = args.shell
    pdb = args.pdb
    outdir = Path(args.d)
    tomap = args.tomap
    res = args.res

    run_particles_mapping(shell=shell, pdb=pdb, outdir=outdir, tomap=tomap, res=res)


if __name__ == "__main__":
    main()