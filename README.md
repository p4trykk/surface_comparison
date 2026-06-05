# SURFMAP
## Requirements
Written in python (3.7 ver) and R (3.6 ver). Algorithm for extraction molecule surface based on MSMS which is already implemented. After modyfication software does not use APBS for electrostatics calculations (instead prepare custom pdb file as input based on PQR file).
Strongly recommended that you install SURFMAP in isolated environment (below short instruction how to create virtual env in different ways).
### How to create isolated environment?
By using isolated virtual environment you will avoid potential version conflicts between python libraries when working on different projects. Below two possible ways how to create and activate virtual env:
1. Create and activate virtual environment (https://pypi.org/project/virtualenv/):
   
   ```
   # upgrade pip to its latest version
   python3 -m pip install --upgrade pip

   # install virtualenv
   python3 -m pip install virtualenv

   # create an isolated environment named 'surfmap_env'
   virtualenv surfmap_env

   # activate isolated environment
   source surfmap_env/bin/activate
   ```
2. Create and activate virtual environment with Conda:
   ```
   # upgrade pip to its latest version
   python3 -m pip install --upgrade pip

   # if you dont have conda - download installer (https://www.anaconda.com/docs/getting-started/miniconda/main) and run
   bash <conda-installer-name>-latest-Linux-x86_64.sh

   # create an isolated environment named 'surfmap_env'
   conda create -n surfmap_env python=3.11

   # activate isolated environment
   conda activate surfmap_env
   ```
Once you're done working on your project, simply type deactivate to exit the environment.

## How to install SURFMAP?
Best way to install SURFMAP is to download project directly from repository. Follow instructions described below if you want modifying the source code:
   ```
   #clone SURFMAP on your machine (type command below in your terminal)
   git clone https://github.com/i2bc/SURFMAP.git

   #go to SURFMAP directory
   cd SURFMAP

   # install SURFMAP
   python3 -m pip install -e .
   ```
### PRoblem with freesasa
If you have problem with FreeSASA during installation, try installing it manually separately: 
```
pip install -U pip setuptools wheel
```
and then the entire project again:
```
pip install -e .
```
If the problem continues, check that version freesas in `setup.py` is correct:
```
freesasa==2.2.1
```
Also make sure that your Python version is min. 3.7 (required, but prefer 3.11).

## How to use
SURFMAP accepts as input a PDB file. Below list of surfmap options:
```
usage: surfmap [-h] (-pdb PDB | -mat MAT | -v) -tomap TOMAP [-proj PROJ] [-res RES] [-rad RAD] [-d D] [-s S] [--nosmooth] [--png] [--keep]
               [--docker] [--pqr PQR] [-ff FF] [-verbose VERBOSE]

options:
  -h, --help            show this help message and exit
  -pdb PDB              Path to a PDB file
  -mat MAT              Input matrix. If the user gives an imput matrix, SURFMAP will directly compute a map from it.
  -v, --version         Print the current version of SURFMAP.
  -tomap TOMAP          Specific key of the feature to map. One of the following: stickiness, kyte_doolittle, wimley_white, electrostatics (not work correctly),
                        circular_variance, bfactor (in our case using as electrostatics), binding_sites, all.
  -proj PROJ            Choice of the projection. Argument must be one of the following: flamsteed, mollweide, lambert. Defaults to flamsteed.
  -res RES              File containing a list of residues to map on the projection. Expected format has the following space/tab separated column
                        values: chainid resid resname
  -rad RAD              Radius in Angstrom added to usual atomic radius (used for calculation solvent excluded surface). The higher the radius the
                        smoother the surface. Defaults to 3.0
  -d D                  Output directory where all files will be written. Defaults to './output_SURFMAP_$pdb_$tomap' with $pdb and $tomap based on
                        -pdb and -tomap given values
  -s S                  Value defining the size of a grid cell. The value must be a multiple of 180. Defaults to 5.0.
  --elec-max-value ELEC_MAX_VALUE
                        Maximum value to be used for the electrostatics color scale. The value will be converted as an absolute value to make the scale symetric around 0. For
                        instance, a value of 5.63 will scale the electrosctatics color values from -5.63 to 5.63.
  --bfactor-min-value BFACTOR_MIN_VALUE
                        Minimum value to be used for the bfactor color scale.
  --bfactor-max-value BFACTOR_MAX_VALUE
                        Maximum value to be used for the bfactor color scale.
  --bfactor-normalize   To normalize scale range to [-1;1].
  --nosmooth            If chosen, the resulted maps are not smoothed (careful: this option should be used only for discrete values!)
  --png                 If chosen, a map in png format is computed (default: only pdf format is generated)
  --keep                If chosen, all intermediary files are kept in the output (default: only final text matrix and pdf map are kept)
  --docker              If chosen, SURFMAP will be run on a docker container (requires docker installed).
  --pqr PQR             Path to a PQR file used for electrostatics calculation. Option only available if '-tomap electrosatics' is requested.
                        Defaults to None.
  -ff FF                Force-field used by pdb2pqr for electrostatics calculation. One of the following: AMBER, CHARMM, PARSE, TYL06, PEOEPB,
                        SWANSON. Defaults to CHARMM.
  -verbose VERBOSE      Verbose level of the console log. 0 for silence, 1 for info level, 2 for debug level. Defaults to 1.
```
You can use pre-made script to automatically launch SURFMAP with predefined settings (`surfmap.sh`). For more information you can go directly to original [SURFMAP repo](https://github.com/i2bc/SURFMAP) or to [SURFMAP publication](https://pubs.acs.org/doi/10.1021/acs.jcim.1c01269).

## General workflow
SURFMAP converts 3D structure from a PDB file into 2D surface map. The standard pipeline begins by generating surface shell around the structure, then assigns value for selected property to each point on that surface. Finally projects these points onto 2D grid and creates smoothed map. 

### PDB as input
The input data is standard PDB file with atom coordinates and B-factor values. For bfactor mode, code uses the x,y,z columns to find atom coordinates and the PDB column corresponding to line[60:67] to find B-factor value. If B-factor value cannot be read, its set to 0. 
In practice, this works as follows: user prepares PDB file in which bfactor column may contain other values for example charges, which allows for efficient creation of electrostatic maps without use of  APBS algorithm.

### Surface creation
SURFMAP does not map all atoms from PDB directly. It first create surface of the structure. The PDB file is converted to the xyzr format, which includes coordinates of atom and radius. Next, MSMS is launched which creates molecular surface (`.vert` files). So map is not map of all atoms, but representation of values assigned to points on the surface.
You can imagine it like this:
```
PDB 3D atoms  →  atoms as spheres  →  molecular surface  →  surface points
```

### 3D to 2D projection
When we assign value to point on the surface, we get spherical coordinates. Then coordinates are projected into surface using coordinate transformations, which are converted to degrees and rounded. The implementation of transformation is different for various types of projections.

### Creating 2D matrix
SURFMAP divides the 2D space into grid cells. The resolution depends on the `-s` parameter which means cell size. Points projected onto the same cells are grouped, and their values are aggregated. Software calculates average value in the cell. Result of this step is matrix in which each cell represents part of the surface and has average value assigned to it.

###Smoothing map
Software smooths values in each cell by averaging values of the nearest neighbors. 



