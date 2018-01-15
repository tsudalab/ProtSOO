#!/bin/bash

# DEFINE ProtPOS HOME DIRECTORY WHERE IT IS INSTALLED
#export PROTPOSHOME="$HOME/opt/protpos/"
export PROTPOSHOME=$HOME/experiment-soo/opt/

# DEFINE GROMACS BINARY PATH 
export GMXBIN="/home/z-hou/soft/gromacs/gmx-double-MPI/bin"
#export GMXBIN="$HOME/opt/gromacs-5.0.7/bin"

# DEFINE PYTHON INTERPRETER
export PYTHONI="/home/z-hou/anaconda24/bin/python"

# DEFINE WHERE TO DO ENERGY MINIMIZATION by GROMACS
export emdirectory="EM"

# DEFINE MOLECULE NAMES 
# FOR PROTEIN, IT CAN BE JUST Protein AS LONG AS GROMACS CAN RECOGNIZE 
# FOR SURFACE, IT SHOULD BE THE RESIDUE NAME 
export proteinm="Protein"
export surfacem="PTS"

# DEFINE INPUT STRUCTURES (PDB)
export protein="protein_lyz.pdb"        
export surface="surface_only.pdb"       

# DEFINE SIMULATION BOX SIZE FOR GROMACS EM, PLEASE ONLY GIVES THREE NUMBERS
# FOR X, Y, and Z DIMENSION IN THE UNIT OF NANOMETER (NM)
export sysboxs="14 14 14"

