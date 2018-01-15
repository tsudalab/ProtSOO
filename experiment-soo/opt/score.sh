#!/bin/bash

#INPUT: $1
#SAMPLE: /absolute/path/to/testcase/pso_conf/conf00001.pdb

echo "================================================================================"
echo "=                                  SCORE.SH                                    ="
echo "================================================================================"

[ -f "set-up.sh" ] && source "set-up.sh" || { echo "Can not find set-up.sh! Please check it!"; exit 1; }

echo "Current working path: "`pwd`
[ -d $emdirectory ] && cd ./$emdirectory || { echo "$emdirectory does not exist!! Maybe wrong input??"; exit 1; }
echo "cd "`pwd`

if [ ! "$1" ]
then echo "Please give the pdb file name!" && exit 1
fi
if [ `echo $1 | cut -c 1` != '/' ]
then echo "the pdb file name should be realpath!" && exit 1
fi


[ -f $PROTPOSHOME/functions.sh ] && source $PROTPOSHOME/functions.sh || { echo "Something wrong with functions.sh! Please check it!"; exit 1; }
checkGromacs    #from functions.sh


export GMX_MAXBACKUP="-1"
CMD=$GMX"editconf"
$GMXBIN/$CMD -box $sysboxs -f "$1" -o "$1"
rm -vf traj.trr traj.xtc ener.edr energy.xvg topol.tpr mdout.mdp md.log confout.gro
CMD=$GMX"grompp"
$GMXBIN/$CMD -c "$1" -f em.mdp -p topol.top -o -maxwarn 5
[ "$?" -ne 0 ] && echo "[grompp] RETURN VALUE: $?" && exit
CMD=$GMX"mdrun"
#$GMXBIN/$CMD -v -s -pd
mpijob -mpi 72  -smp 1  $GMXBIN/$CMD -v -s
#$GMXBIN/$CMD -v -s 
[ "$?" -ne 0 ] && echo "[mdrun] RETURN VALUE: $?" && exit
GENG=g_energy
if [ ! -x $GMXBIN/$GENG ]
then
    GENG="gmx energy"
fi
#ASSUMED the name is first protein then following by surface
courn=`echo "" | $GMXBIN/$GENG -f 2>&1 | grep -Po "([0-9]+)(?=\s+Coul-SR:$proteinm-$surfacem)" | head -n 1`
ljn=`echo "" | $GMXBIN/$GENG -f 2>&1 | grep -Po "([0-9]+)(?=\s+LJ-SR:$proteinm-$surfacem)" | head -n 1`
echo "We will select energy number: $courn and $ljn"
echo "Please confirm they are correct."
echo -e "$courn\n$ljn" | $GMXBIN/$GENG -f ener.edr -o energy.xvg

echo "================================================================================"

#OUTPUT file: energy.xvg (under EM directory)
#SAMPLE:
#   54.000000  -197.174576  -49.886299
