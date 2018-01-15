#!/bin/bash


function checkInstall   {
    checkPython

    checkGromacs

    #check grep perl expression supporting
    { echo "perl expression testing" | grep -Po '(?=.*)'; } || { echo "Please install GNU grep!"; exit 1; }

    echo "Installation checking passed."
}

function checkPython    {
    for mod in "pymol" "numpy" "sklearn"
    do
        echo "import $mod" | $PYTHONI || { echo "Please install $mod"; exit 1; }
    done
}

function checkGromacs   {
    echo "Gromacs PATH: $GMXBIN"
    GMX=""

    if [ -x $GMXBIN/gmx ]
    then
        GMX="gmx "
    elif [ -x $GMXBIN/editconf ] && [ -x $GMXBIN/grompp ] && [ -x $GMXBIN/mdrun ]
    then
        echo "Gromacs environment confirmed."
    elif [ -x /usr/bin/editconf ] && [ -x /usr/bin/grompp ] && [ -x /usr/bin/mdrun ]
    then
        echo "Using the system installed gromacs."
        GMXBIN="/usr/bin"
    else
        echo "Gromacs environment invalid!!"
        echo "EXIT NOW!!"
        exit 1
    fi
}

function checkenv {
    #ASSUMED python is corrent, current version will not check that.
    chklst=("$PROTPOSHOME/simplePSO.py" "$PROTPOSHOME/simpleANS.py" "$PROTPOSHOME/simpleMOVE.py" $PYTHONI)
    for i in ${chklst[*]}
    do
        [ -f $i ] || { echo "Can not find $i!!"; exit 1; }
    done

    gmxv=4
    [ -x $GMXBIN/gmx ] && gmxv=5
    #It should do only once for each testing
    #echo $i | rev | cut -d '.' -f 2- | rev
    pred=`pwd`
    cd $emdirectory || { echo "emdirectory: $emdirectory does not exist! Please check set-up.sh!"; exit 1; }
    echo "prepare the mdp file from template"
    if [ -f "em.mdp" ]
    then
        echo "em.mdp already exist. Created a backup"
        datecode=`date '+%Y%m%d%H%M%S%N'`
        mv -v em.mdp em.mdp.$datecode
    fi
    awk  -v pname="$proteinm" -v sname="$surfacem" -v gmxv=$gmxv \
    '
    /^\s*energygrps\s+=/{print "energygrps = "pname"  "sname; next;}
    /^\s*energygrp_excl\s+=/{print "energygrp_excl = "sname"  "sname; next;}
    /^\s*freezegrps\s+=/{print "freezegrps = "sname; next;}
    /^\s*cutoff-scheme\s+=/{if (gmxv=5) { print $0; }; next;}
    1
    ' em.mdp.tpl > em.mdp
    emf=em.mdp.$datecode
    echo "Duplicated backup will be removed."
    diff -q em.mdp $emf && rm -vf $emf  #remove if the contents are the same
    cd $pred
}


function simplepso {
    checkenv

    flist=({debug,score}.log db.json   gbest{.pdb,.txt,_sorted.txt,_energy.txt,_vector.txt} )
    echo "=============================================================="
    echo "    ProtPOS STARTED @ "`date '+%Y-%m-%d  %H:%M:%S  '`
    echo "=============================================================="
    echo "removing the previous run output files "
    rm -rf ${flist[*]}
    echo "the protein is: $protein"
    echo "the surface is: $surface"
    echo "=============================================================="
    $PYTHONI $PROTPOSHOME/te.py "$@"  
    [ $? -ne 0 ] && exit 1    #stop if the core searching procedure is failed
    ####### ANALYSISING RESULTS #############################################################################################################
    $PYTHONI $PROTPOSHOME/simpleANS.py --pdb gbest.pdb --jdb db.json
    [ $? -ne 0 ] && exit 1    #stop if the core searching procedure is failed
    echo "=============================================================="
    datn=`date '+%m%d%H%M'`
    [ -e "protpos-$datn" ] && echo "WARNING!! protpos-$datn is already existed! Files may be overwritten!"
    mkdir "protpos-$datn"
    mv    ${flist[*]} "protpos-$datn"
    [ "$?" -ne 0 ] && echo "Can not copy some files, the output may be corrupted"
    echo "Packed the run result data into directory: protpos-$datn"
    echo "=============================================================="
    echo "    "`date '+%Y-%m-%d  %H:%M:%S  '`" @ ProtPOS END    "
    echo "=============================================================="
}

function clustering {
    DBCLUSTER=$PROTPOSHOME/dbcluster.py
    [ ! -x $DBCLUSTER ] && { echo "Please check $PAOSHOME!"; exit 1; }

    [ $1   == "" ] &&  $1=6.0
    [ -x ./cluster ] && mv -v cluster cluster-`date '+%Y%m%d%H%M%S%N'`
    find -L ./ -iname 'gbest.txt'  -exec $PYTHONI $DBCLUSTER cluster  --out-dir "cluster"  --minn 2 --eps $1   --gbs {} +

    mv -v debug.log cluster/cluster.log
    cn=`ls cluster/cluster*.pdf | wc -l`
    echo "$cn cluster(s) were found."
    echo "For more details, please look at cluster directory."
}
