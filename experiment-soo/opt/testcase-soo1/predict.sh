#!/bin/bash

cdr=`pwd`
echo "Current working directory is: $cdr"


# SETTING UP THE ENVIRONMENT
[ -f $cdr/set-up.sh ] && source $cdr/set-up.sh || { echo "Can not find set-up.sh in the current working directory!"; exit 1; }
[ -d "$PROTPOSHOME" ] || { echo "$PROTPOSHOME does not exist!!"; exit 1; }
[ -f $PROTPOSHOME/functions.sh ] && source $PROTPOSHOME/functions.sh || { echo "Something wrong with functions.sh! Please check it!"; exit 1; }


# VERIFYTING PROTPOS INSTALLATION AND REQUIRED SOFTWARE
checkInstall    


# REQUIRED PARAMETERS FOR SIMPLEPSO
# --maxx (angstrom) Upper limit for translation in X direction 
# --minx (angstrom) Lower limit for translation in X direction 
# --maxy (angstrom) Upper limit for translation in Y direction 
# --miny (angstrom) Lower limit for translation in Y direction 

# OPTIONAL PARAMETERS
# --r    Convergence criteria (default=10)
# --n    Number of particles  (default=200)
# --c1   Congitive weight     (default=1.193)
# --c2   Social weight        (default=1.193)
# --w    Inertia weight       (default=0.721)
# --maxz (angstrom) Upper limit for translation in Z direction relative to the surface (default=5.5)
# --minz (angstrom) Lower limit for translation in Z direction relative to the surface (default=1.0)

#--resi   Find protein orientations containing any of the specified contacting residues.
#         e.g. residue ID 10 or 20: --resi 10 20
#              residue ID 10 to 15: --resi {10..15}

#--init   If set, protein position and orientation with respect to the surface are used as the initial
#         structure for the search (default is unset, means position at center of surface and random 
#         orientation). This feature helps to force sampling specific region of the surface

#--offset (decimal, in format Rx Ry Rz Tx Ty Tz) generate the initial structure by translating and 
#         rotating the given protein structure instead of a random orientation


# RUN PSO TO PREDICT PREFERRED PROTEIN ORIENTATION ON SURFACE

#==================================================================================
# RUN ONCE -- uncomment the following 2 lines
#==================================================================================
 simplepso --proteinf=$protein --surfacef=$surface --r 10 --n 2 --c1 1.193  --c2 1.193  --w 0.721 \
          


#==================================================================================
# RUN 10 TIMES AND PERFORM CLUSTERING -- uncomment the following 8 lines 
# NOTE WE USED A SMALLER NUMBER OF PSO PARTICLES (n) SO RUNTIME WILL BE SHORTER
# BUT LESS FAVORABLE ORIENTATIONS MAY BE FOUND COMPARED TO PUBLISHED RESULTS
#==================================================================================
# for i in {1..10};
# do
# simplepso --proteinf=$protein --surfacef=$surface --r 10 --n 100 --c1 1.193  --c2 1.193  --w 0.721 \
#           --maxx=+3.0 --minx=-3.0 --maxy=+3.0 --miny=-3.0 --maxz=5.5 --minz=1.0 
# done
#
# EPS=6.0         
# clustering $EPS


#==================================================================================
# LOOKING FOR ORIENTATION WITH SPECIFIC CONTACTING RESIDUES
#                                                -- uncomment the following 2 lines
#===================================================================================
# simplepso --proteinf=$protein --surfacef=$surface --r 10 --n 200 --c1 1.193  --c2 1.193  --w 0.721 \
#           --maxx=+3.0 --minx=-3.0 --maxy=+3.0 --miny=-3.0 --maxz=5.5 --minz=1.0 --resi 14 15   
 

#==================================================================================
# LOOKING FOR ORIENTATION AT CERTAIN REGION OF SURFACE  
#                                                -- uncomment the following 3 lines
# NOTE WE USED A DIFFERENT PROTEIN PDB HERE IN WHICH THE PROTEIN IS AT A SPECIFIC
# POSITION OF THE SURFACE
#===================================================================================
# export protein="protein_lyz_atcorner.pdb"
# simplepso --proteinf=$protein --surfacef=$surface --r 10 --n 200 --c1 1.193  --c2 1.193  --w 0.721 \
#           --maxx=+3.0 --minx=-3.0 --maxy=+3.0 --miny=-3.0 --maxz=5.5 --minz=1.0 --init 


