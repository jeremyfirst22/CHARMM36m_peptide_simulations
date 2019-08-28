#!/bin/bash
## Jeremy First

TOP=${PWD}
MDP=$TOP/mdp_files
FF=$TOP/GMXFF
forceField=charmm36
saltConcentration=0
numMols=14
spacing=0.497 ##nm 
glyDist=0.6822 ##nm from geometry optimization, Spartan '16, DFT wB97X-D, 6-311+g**
glyRest=1000   #kJ/mol/nm
sulRest=200000 #kJ/mol/nm
fileName=$TOP/StartingStructures/folded.pdb
totSimTime=50
zBuff=0.5 #nm. This is the spacing of the SAM layer off the z intercept. Must not be zero, 
          ##    or the SAM splits during NPT box growth. 

verbose=false
analysis=false

usage(){
    echo "USAGE: $0 -f prep/folded/unfolded -L C11/LIG [ options ]" 
    exit 
}

HELP(){
    echo 
    echo "This program runs molecular dynamic simulations of a peptide on a SAM surface" 
    echo 
    echo "Usage: $0 [options] "
    echo "  -f   folded, unfolded, or prep. Mandatory" 
    echo "  -L   LIG or C11. Mandtory. Only LIG/C11 SAM molecules supported currently." 
    echo "  -c   Starting structure  (folded) PDB file: Default = $fileName "
    echo "  -t   Maximum simulation time (ns) : Default = $totSimTime "
    echo "  -a   Perform analyis on trajectory : Default = no"
    echo "  -s   Salt conentration (mM) : Default = $saltConcentration (neutralizing)" 
    echo "  -n   Number of decanethiol molecules in each direction of the layer : Default = $numMols" 
    echo "  -D   Distance (nm) between each decanethiol molecule : Default = $spacing" 
    echo "  -g   Distnace (nm) between the alpha carbon of glycines to C10 of the SAM layer" 
    echo "              Default = $glyDist from Geoemetry optimization, Spartan '16, DFT wB97X-D, 6-311+g**"
    echo "  -r   Restraint force constant (kJ/mol/nm) of the "bond" from glycine to nearest decanethiol C10"
    echo "              Default = $glyRest "
    echo "  -R   Restraint force constant (kJ/mol/nm) of the "bond" keeping sulfur atoms in plane" 
    echo "              Default = $sulRest " 
    echo "  -z   Distance of the bottom of the SAM layer to the z-axis of the simulation box. " 
    echo "              Default = $zBuff  " 
    echo "  -m   Location of the mdp_files : Default = $MDP"
    echo "  -p   Location of force field files. : Default = $FF"
    echo "  -n   Name of force field : Default = $forceField"
    echo "  -v   Print all options and quit." 
    echo "  -h   print this usage and exit "
    echo "" 
    exit
}

while getopts :f:L:c:t:as:n:D:g:r:R:z:m:p:n:vh opt; do 
   case $opt in 
      f) 
        fold=$OPTARG
        ;; 
      L) 
        RESNAME=$OPTARG
        ;; 
      c)
        fileName=${TOP}/$OPTARG
        ;; 
      t)
        totSimTime=$OPTARG
        ;; 
      a) 
        analysis=true
        ;; 
      s) 
        saltConcentration=$OPTARG
        ;; 
      D) 
        spacing=$OPTARG
        ;; 
      g) 
        glyDist=$OPTARG
        ;; 
      r) 
        glyRest=$OPTARG
        ;; 
      R) 
        sulRest=$OPTARG
        ;; 
      z) 
        zBuff=$OPTARG
        ;; 
      m) 
        MDP=${TOP}/$OPTARG
        ;; 
      p) 
        FF=${TOP}/$OPTARG
        ;; 
      n) 
        forceField=$OPTARG
        ;; 
      v) 
        verbose=true
        ;; 
      :) 
        echo " option -$OPTARG requires an argument! "
        usage
        ;; 
      h)
        HELP
        ;; 
      \?) 
        echo "Invalid option -$OPTARG "
        HELP 
        ;; 
   esac
   done 

main(){
    checkInput
    printf "\n\t\t*** Program Beginning $SOL_$fold $totSimTime (ns)***\n\n" 
    if [ ! -d $SOL ] ; then mkdir $SOL ; fi 
    cd $SOL
    if [ $fold == "prep" ] ; then 
        prep
    else  
        production 
        if $analysis ; then 
            analysis
            fi 
        fi 
    cd ../
    printf "\n\n\t\t*** Program Ending    ***\n\n" 
}

prep(){
    if [ ! -d prep ] ; then mkdir prep ; fi 
    cd prep
    build_SAM
    layer_relax
    protein_steep
    solvate
    solvent_steep
    solvent_nvt
    solvent_npt
    build_system
    system_steep
    system_nvt
    system_npt
    heating
    heated
    cooling
    cd ../
}

analysis(){
    if [ ! -d $fold ] ; then mkdir $fold ; fi 
    cd $fold 
    #dssp
    #rgyr
    #minimage
    #rdf
    #nopbc
    #rmsd 
    #cd_spectra
    #cluster
    #good-turing
    #order
    density 
    hbond
    cd ../
}

checkInput(){
    if [ $saltConcentration -eq 0 ] ; then 
        SOL=${RESNAME}_sam
    else 
        SOL=${RESNAME}_sam_$saltConcentration
    fi 
    MOLEC=folded_$SOL
    logFile=$TOP/$SOL/$fold/$fold.log
    errFile=$TOP/$SOL/$fold/$fold.err
    if $verbose ; then 
        echo "MOLEC : $MOLEC"
        echo "RESNAME : $RESNAME"
        echo "SOL : $SOL"
        echo "Folded state : $fold" 
        echo "Input file name: $fileName"
        echo "Max simultaiton time: $totSimTime"
        echo "Salt concentration : $saltConcentration"
        echo "Number of decanethiols: $numMols" 
        echo "Distance between each decanethiol : $spacing"
        echo "Distance of glycine to SAM : $glyDist" 
        echo "Glycine restraint : $glyRest"
        echo "Sulfur restraint : $sulRest" 
        echo "Perform analysis : $analysis " 
        echo "Verbose = $verbose"
        echo "Distance to z-axis: $zBuff" 
        echo "mdp files = $MDP " 
        echo "force field directory = $FF" 
        echo "force field name = $forceField " 
        echo "Log file = $logFile " 
        echo "Error file = $errFile " 
        echo ""
        exit 
        fi 
    
    if [ ! -f $fileName ] ; then 
        echo "ERROR: $fileName not found " 
        exit 
        fi 
    if [[ $fileName != *.pdb ]] ; then 
        echo "ERROR: Input file must be PDB file (*.pdb)" 
        exit 
        fi 
    if [ ! -z $fold ] ; then 
        if [[ $fold != 'folded' && $fold != 'unfolded' && $fold != 'prep' ]] ; then 
            echo "ERROR: $fold not recognized"
            HELP
            fi 
    else 
        echo "ERROR: Must specify folded/unfoled" 
        usage
        fi 
    if [ ! -z $RESNAME ] ; then 
        if [[ $RESNAME == 'C11' ]] ; then 
            RESFIRSTH=HS
            RESSUL=S1
            RESLASTHEAVY=O24
        elif [[ $RESNAME == 'LIG' ]] ; then 
            RESFIRSTH=H22
            RESSUL=S1
            RESLASTHEAVY=C10
        else 
            echo "ERROR: $RESNAME not recognized"
            HELP
            fi 
        decanethiolFile=$TOP/StartingStructures/$RESNAME.pdb
        if [ ! -f $decanethiolFile ] ; then 
            echo "ERROR: $decanethiolFile not found " 
            exit 
            fi 
    else 
        echo "ERROR: Must specify SAM molecule name: C11/LIG" 
        usage
        fi 
    if [ $((numMols%2)) -ne 0 ] ; then 
        echo "ERROR: numMols must be even for sam layer to tesselate" 
        exit 
        fi 
    if [[ $glyDist < 0 ]] ; then 
        echo "ERROR: Glycine distance to SAM must be larger than zero" 
        exit 
        fi 
    if [[ $spacing < 0 ]] ; then 
        echo "ERROR: spacing between dcanethiols must be larger than zero" 
        exit 
        fi 
    if [[ $glyRest < 0 ]] ; then 
        echo "ERROR: glycine restraint force constant must be larger than zero" 
        exit 
        fi 
    if [[ $sulRest < 0 ]] ; then 
        echo "ERROR: sulfur restraint force constant must be larger than zero" 
        exit 
        fi 
    if [ $totSimTime -le 0 ] ; then 
        echo "ERROR: Total simulation time must be larger than 0" 
        exit 
        fi 
    if [ ! -d $FF/$forceField.ff ] ; then 
        echo ; echo "ERROR: FF not found" 
        exit
        fi  
    if [ ! -d $MDP ] ; then 
        echo ; echo "ERROR: mdp files directory not found" 
        exit 
        fi 
    if [[ $saltConcentration < 0 ]] ; then 
        echo "ERROR: Negative salt conentration given" 
        exit 
        fi 
    check $MDP $fileName $FF $FF/$forceField.ff 
} 

check(){
   for arg in $@ ; do 
        if [ ! -s $arg ] ; then 
            echo ; echo "ERROR: $arg missing. Exitting" 
            exit 
            fi 
        done 
}

clean(){
    if [ -d $forceField.ff ] ; then rm -r $forceField.ff *.dat ; fi 
}

create_dir(){
    if [ -z $1 ] ; then 
        echo "ERROR: create_dir requires argument. " ; exit ; fi 

    dirName=$1 
    if [ ! -d $dirName ] ; then mkdir $dirName ; fi 
    
    if [ ! -d $dirName/$forceField.ff ] ; then 
        if [ -d $FF/$forceField.ff ] ; then 
            cp -r $FF/$forceField.ff $dirName
            cp $FF/*.dat $dirName/. 
        else 
            echo "FF not found" 
            exit 
            fi 
        fi 
}

build_SAM(){
    printf "\t\tBuilding SAM layer........................" 
    if [ ! -f Build_SAM/bottom.gro ] ; then 
        create_dir Build_SAM 
        
        #cp $MOLEC.pdb Build_SAM/. 
        cp $decanethiolFile Build_SAM/.
        cd Build_SAM

        echo 'System' | gmx editconf -f `basename $decanethiolFile` \
            -princ \
            -o oriented.gro >> $logFile 2>> $errFile 
        check oriented.gro 

        gmx editconf -f oriented.gro \
            -rotate 180 90 0 \
            -o upright.gro >> $logFile 2>> $errFile 
        check upright.gro 

        gmx editconf -f upright.gro \
            -rotate 30 0 0 \
            -center 0 0 0 \
            -o inplane.gro >> $logFile 2>> $errFile  
        check inplane.gro 

        gmx editconf -f inplane.gro \
            -rotate 0 22 0 \
            -center 0 0 0 \
            -o twisted.gro >> $logFile 2>> $errFile 
        check twisted.gro 

        echo "#!/usr/bin/env python
import math as m

numMols=$numMols
spacing=$spacing
step=spacing * m.cos(m.pi / 6 ) 

for i in range(numMols) : 
    for j in range(numMols) : 
        x,y,z = i*step, j*spacing, 0.0 
        if not i%2 == 0 : 
             y += spacing / 2
        print \"%.3f  %.3f  %.3f\"%(x,y,z) 
" > make_position.py 

        python make_position.py > position.dat 
        check position.dat 
        
        xdim=`echo "$numMols * $spacing * c(4*a(1) / 6) " | bc -l | awk '{printf "%f", $0}'`
        ydim=`echo "$numMols * $spacing " | bc -l | awk '{printf "%f", $0}'`
        zdim=$xdim

        gmx insert-molecules -ci twisted.gro \
            -ip position.dat \
            -rot none \
            -box $xdim $ydim $zdim \
            -o layer.gro >> $logFile 2>> $errFile 
        check layer.gro 

        gmx editconf -f layer.gro \
            -bt triclinic \
            -o boxed.gro >> $logFile 2>> $errFile  
        check boxed.gro 

        gmx pdb2gmx -f boxed.gro \
            -p layer.top \
            -ff $forceField \
            -water tip3p \
            -o boxed.gro >> $logFile 2>> $errFile 
        check layer.top 

        includeLine=`cat -n layer.top | grep '; Include Position restraint file' | awk '{print $1}' | tail -n1 `
        head -n $includeLine layer.top > decanethiol.top 
        echo '#ifdef POSSULRES 
#include "posre_SUL.itp" 
#endif' >> decanethiol.top 
        ((includeLine++)) 
        tail -n +$includeLine layer.top >> decanethiol.top

        echo "[ position_restraints ]" > posre_SUL.itp 
        echo ";; Pin sulfur atoms to spacing found on SAM" >> posre_SUL.itp 
        for atom in `grep $RESNAME boxed.gro  | cut -c9- | grep -E "^[ ]*$RESSUL" | cut -c8- | awk '{print $1}'` ; do 
            printf "%6i%6i%10.f%10.f%10.f\n" $atom 1 0 0 $sulRest >> posre_SUL.itp 
        done 

        zshift=`grep $RESNAME boxed.gro | grep " $RESFIRSTH " | awk '{print $6}' | sort -n | uniq | head -n1`
        yshift=`grep $RESNAME boxed.gro | grep " $RESFIRSTH " | awk '{print $5}' | sort -n | uniq | head -n1`
        xshift=`grep $RESNAME boxed.gro | grep " $RESFIRSTH " | awk '{print $4}' | sort -n | uniq | head -n1`
    
        zshift=`echo "-1 * $zshift " | bc -l | awk '{printf "%f", $0}'`
        yshift=`echo "-1 * $yshift " | bc -l | awk '{printf "%f", $0}'`
        xshift=`echo "-1 * $xshift " | bc -l | awk '{printf "%f", $0}'`

        gmx editconf -f boxed.gro \
            -translate $xshift $yshift $zshift \
            -o bottom.gro >> $logFile 2>> $errFile 
        check bottom.gro

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
    check Build_SAM/posre_SUL.itp Build_SAM/decanethiol.top Build_SAM/bottom.gro 
} 

layer_relax(){
    printf "\t\tSteep and NVT relax of SAM layer.........." 
    if [ ! -f Relax_SAM/nvt_relax.nopbc.gro ] ; then 
        create_dir Relax_SAM
        
        cp Build_SAM/decanethiol.top Relax_SAM/.
        cp Build_SAM/bottom.gro Relax_SAM/.
        cp Build_SAM/*.itp Relax_SAM/. 
        cd Relax_SAM
        
        gmx grompp -f $MDP/sam_steep.mdp \
            -c bottom.gro \
            -p decanethiol.top \
            -o steep.tpr >> $logFile 2>> $errFile 
        check steep.tpr 

        gmx mdrun -deffnm steep >> $logFile 2>> $errFile 
        check steep.gro 

        gmx grompp -f $MDP/sam_nvt_relax.mdp \
            -c steep.gro \
            -p decanethiol.top \
            -o nvt_relax.tpr >> $logFile 2>> $errFile 
        check nvt_relax.tpr 

        if [ ! -f nvt_relax.gro ] ; then 
            gmx mdrun -deffnm nvt_relax >> $logFile 2>> $errFile 
            fi 
        check nvt_relax.gro 

        echo 'System' | gmx trjconv -f nvt_relax.gro \
            -s nvt_relax.tpr \
            -pbc nojump \
            -o nvt_relax.nopbc.gro >> $logFile 2>> $errFile 
        check nvt_relax.nopbc.gro

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

protein_steep(){
    printf "\t\tProtein steep............................." 
    if [ ! -f Protein_steep/protein_steep.gro ] ; then 
        create_dir Protein_steep
        
        cp $fileName Protein_steep/.
        cd Protein_steep

        ## 3, 4 -- None, None for termini options
        echo '3 4' | gmx pdb2gmx -f $(basename $fileName) \
            -p $MOLEC.top \
            -ff $forceField \
            -ter \
            -water tip3p \
            -o $MOLEC.gro >> $logFile 2>> $errFile 
        check $MOLEC.gro 

        xdim=`tail -n1 ../Relax_SAM/nvt_relax.gro | awk '{print $1}'`
        ydim=`tail -n1 ../Relax_SAM/nvt_relax.gro | awk '{print $2}'`
        zdim=`tail -n1 ../Relax_SAM/nvt_relax.gro | awk '{print $3}'`

        echo 'Backbone' | gmx editconf -f $MOLEC.gro \
            -box $xdim $ydim $zdim \
            -bt tric \
            -princ \
            -center 0 0 0 \
            -o aligned.gro >> $logFile 2>> $errFile
        check aligned.gro 

        ## 90 deg around z axis puts long axis in long dimension of box
        ## 0 deg around x axis puts CA of Gly at bottom, with Leus closest to SAM
        gmx editconf -f aligned.gro \
            -rotate 0 0 90 \
            -o rotated.gro >> $logFile 2>> $errFile 
        check rotated.gro 

        ##This is the distance to the edge of the box. 
        xshift=`echo "$xdim / 2" | bc -l | awk '{printf "%f", $0}'`
        yshift=`echo "$ydim / 2" | bc -l | awk '{printf "%f", $0}'`

        zshift=`grep GLY rotated.gro | grep CA | awk '{total += $6} END {print total/NR}'`
        zshift=`echo "$zshift * -1 + $glyDist - 0.12" | bc -l | awk '{printf "%f", $0}'`
        gmx editconf -f rotated.gro \
            -translate $xshift $yshift $zshift \
            -o translated.gro >> $logFile 2>> $errFile  
        check translated.gro 

        gmx grompp -f $MDP/protein_steep.mdp \
            -c translated.gro \
            -p $MOLEC.top \
            -o protein_steep.tpr >> $logFile 2>> $errFile 
        check protein_steep.tpr 

        gmx mdrun -deffnm protein_steep >> $logFile 2>> $errFile 
        check protein_steep.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi
} 

solvate(){
    printf "\t\tSolvating protein........................." 
    if [ ! -f Solvate/neutral.top ] ; then 
        create_dir Solvate
        
        cp Protein_steep/protein_steep.gro Solvate/. 
        cp Protein_steep/$MOLEC.top Solvate/. 
        cd Solvate

        xdim=`tail -n1 protein_steep.gro | awk '{print $1}'`
        ydim=`tail -n1 protein_steep.gro | awk '{print $2}'`
        zdim=`tail -n1 protein_steep.gro | awk '{print $3}'`

        gmx solvate -cp protein_steep.gro \
            -box $xdim $ydim $zdim \
            -p $MOLEC.top \
            -o solvated.gro >> $logFile 2>> $errFile 
        check solvated.gro

        gmx grompp -f $MDP/vac_md.mdp \
            -p $MOLEC.top \
            -c solvated.gro \
            -o genion.tpr >> $logFile 2>> $errFile 
        check genion.tpr
        
        protCharge=$(grep qtot $MOLEC.top | tail -n -1 | awk '{print $11}') 
        if [[ $protCharge > 0 ]] ; then 
            numNa=$(python -c "print int(round($saltConcentration * 6.022*10**-4 * $xdim*$ydim*$zdim))") 
            numCl=$(python -c "print $protCharge + $numNa") ##neutral box required
        else   
            numCl=$(python -c "print int(round($saltConcentration * 6.022*10**-4 * $xdim*$ydim*$zdim))") 
            numNa=$(python -c "print -1 * $protCharge + $numCl") ##neutral box required
        fi 

        echo 'Water' | gmx genion -s genion.tpr \
            -nn $numCl \
            -np $numNa \
            -nname 'CLA' \
            -pname 'SOD' \
            -o neutral.gro >> $logFile 2>> $errFile 
        check neutral.gro 

        ## 3, 4 -- None, None for terimini options
        echo '3 4' | gmx pdb2gmx -f neutral.gro \
            -ff $forceField \
            -water tip3p \
            -p neutral.top \
            -ter \
            -o neutral.gro >> $logFile 2>> $errFile 
        check neutral.top 

        sed 's/POSRES/POSRES_IONS/' neutral_Ion2.itp > temp.itp 
        mv temp.itp neutral_Ion2.itp 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

solvent_steep(){
    printf "\t\tSolvent steep............................." 
    if [ ! -f Solvent_steep/solvent_steep.gro ] ; then 
        create_dir Solvent_steep
        
        cp Solvate/neutral.gro Solvent_steep/. 
        cp Solvate/neutral.top Solvent_steep/. 
        cp Solvate/neutral_*.itp Solvent_steep/. 
        cp Solvate/posre_*.itp Solvent_steep/. 
        cd Solvent_steep

        gmx grompp -f $MDP/solvent_steep_easy.mdp \
            -p neutral.top \
            -c neutral.gro \
            -o solvent_steep_easy.tpr >> $logFile 2>> $errFile 
        check solvent_steep_easy.tpr 

        gmx mdrun -deffnm solvent_steep_easy -nt 1 >> $logFile 2>> $errFile 
        check solvent_steep_easy.gro 

        gmx grompp -f $MDP/solvent_steep.mdp \
            -p neutral.top \
            -c solvent_steep_easy.gro \
            -o solvent_steep.tpr >> $logFile 2>> $errFile 
        check solvent_steep.tpr 

        gmx mdrun -deffnm solvent_steep >> $logFile 2>> $errFile 
        check solvent_steep.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

solvent_nvt(){
    printf "\t\tSolvent NVT relaxation...................." 
    if [ ! -f Solvent_nvt/solvent_nvt.gro ] ; then 
        create_dir Solvent_nvt
        
        cp Solvent_steep/solvent_steep.gro Solvent_nvt/. 
        cp Solvent_steep/neutral.top Solvent_nvt/. 
        cp Solvent_steep/*.itp Solvent_nvt/. 
        cd Solvent_nvt

        gmx grompp -f $MDP/solvent_nvt_relax.mdp \
            -c solvent_steep.gro \
            -p neutral.top \
            -o solvent_nvt.tpr >> $logFile 2>> $errFile 
        check solvent_nvt.tpr 

        gmx mdrun -deffnm solvent_nvt >> $logFile 2>> $errFile 
        check solvent_nvt.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

solvent_npt(){
    printf "\t\tSolvent NPT semiisotropic relaxation......" 
    if [ ! -f Solvent_npt/solvent_npt.gro ] ; then 
        create_dir Solvent_npt
        
        cp Solvent_nvt/solvent_nvt.gro Solvent_npt/. 
        cp Solvent_nvt/neutral.top Solvent_npt/. 
        cp Solvent_nvt/*.itp Solvent_npt/. 
        cd Solvent_npt

        gmx grompp -f $MDP/solvent_semiisotropic_npt_relax.mdp \
            -c solvent_nvt.gro \
            -p neutral.top \
            -o solvent_npt.tpr >> $logFile 2>> $errFile 
        check solvent_npt.tpr 

        gmx mdrun -deffnm solvent_npt >> $logFile 2>> $errFile 
        check solvent_npt.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

build_system(){
    printf "\t\tBuilding surface+protein+water system....." 
    if [ ! -f Build_system/system.top ] ; then 
        create_dir Build_system 
        
        cp Solvent_npt/solvent_npt.gro Build_system/. 
        #cp Solvent_npt/neutral.top Build_system/. 
        #cp Solvent_npt/*.itp Build_system/. 
        cp Relax_SAM/nvt_relax.nopbc.gro Build_system/. 
        cd Build_system

        zdim=`tail -n1 solvent_npt.gro | awk '{print $3}'`
        ydim=`tail -n1 solvent_npt.gro | awk '{print $2}'`
        xdim=`tail -n1 solvent_npt.gro | awk '{print $1}'`
        zshift=`cat nvt_relax.nopbc.gro | grep $RESNAME | grep $RESLASTHEAVY | awk '{print $6}' | sort -n | tail -n1`
    
        gmx editconf -f nvt_relax.nopbc.gro \
            -translate 0 0 $zBuff \
            -o bottom_boxed.gro >> $logFile 2>> $errFile 
        check bottom_boxed.gro
    
        zdim=`echo "$zdim + $zshift + 0.05" | bc -l | awk '{printf "%f", $0}'`

        gmx editconf -f solvent_npt.gro \
            -box $xdim $ydim $zdim \
            -o new_box.gro >> $logFile 2>> $errFile 
        check new_box.gro 

        zshift=`cat new_box.gro | grep SOL | grep OW | awk '{print $6}' | sort -n | uniq | tail -n1`
        zshift=`echo "$zdim - $zshift + $zBuff + 0.05" | bc -l | awk '{printf "%f", $0}'` ##1 angstrom buffer space to 
                                                                                         ##  keep atoms from being on top of each other
    
        gmx editconf -f new_box.gro \
            -translate 0 0 $zshift \
            -o shifted.gro >> $logFile 2>> $errFile 
        check shifted.gro 

        head -n1 shifted.gro > system.gro 
        protAtoms=`cat shifted.gro | tail -n +2 | head -n 1`
        samAtoms=`cat bottom_boxed.gro | tail -n +2 | head -n 1`
        totAtoms=`echo "$protAtoms + $samAtoms" | bc -l | awk '{printf "%f", $0}'`

        echo $totAtoms >> system.gro 
        tail -n +3 shifted.gro | sed '$ d' | grep -v 'SOL' | grep -v 'CLA' | grep -v 'SOD' >> system.gro 
        tail -n +3 bottom_boxed.gro | sed '$ d' >> system.gro 
        tail -n +3 shifted.gro | grep 'CLA' >> system.gro 
        tail -n +3 shifted.gro | grep 'SOD' >> system.gro 
        tail -n +3 shifted.gro | grep 'SOL' >> system.gro 
        tail -n1 shifted.gro >> system.gro 

        ## 3, 4 -- None, None for terimini options
        echo '3 4' | gmx pdb2gmx -f system.gro \
            -p neutral.top \
            -ter \
            -ff $forceField \
            -water tip3p \
            -merge interactive \
            -o system.gro >> $logFile 2>> $errFile
        check neutral.top 

        includeLine=`cat -n neutral.top | grep '; Include Position restraint file' | awk '{print $1}'`
        head -n $includeLine neutral.top > system.top  
        echo '#include "distance_restraints.itp"' >> system.top 
        echo >> system.top 
        echo '#ifdef POSSULRES' >> system.top 
        echo '#include "posre_SUL.itp"' >> system.top 
        echo '#endif' >> system.top 
        echo >> system.top 
        ((includeLine++)) 
        tail -n +$includeLine neutral.top >> system.top 

        ##Sulfer atoms have changed numbers. Need to rebuild posre_SUL.itp 
        echo "[ position_restraints ]" > posre_SUL.itp 
        echo ";; Pin sulfur atoms to spacing found on SAM" >> posre_SUL.itp 
        for atom in `grep $RESNAME system.gro  | cut -c9- | grep -E "^[ ]*$RESSUL" | cut -c8- | awk '{print $1}'` ; do 
            printf "%6i%6i%10.f%10.f%10.f\n" $atom 1 0 0 $sulRest >> posre_SUL.itp 
        done 

        ## We need to rebuild posre.itp. pdb2gmx assigns all heavy atoms to restraints. 
        ## We want the decanethiol heavy atoms to be able to relax with protein position restraints. 
        echo "; Only protein heavy atom position restraints" > posre.itp 
        echo "[ position_restraints ]" >> posre.itp 
        echo "; atom  type      fx      fy      fz" >> posre.itp 
        for atom in `tail -n +2 system.gro | sed '$ d' | grep -v HOH | grep -v CLA | grep -v SOD | grep -v $RESNAME | awk '{print $2"\t"$3}' | grep "^[CNOS]" | awk '{print $2}'` ; do 
            printf "%6i%6i%10.f%10.f%10.f\n" $atom 1 1000 1000 1000 >> posre.itp 
            done 

        add_restraints

        check system.top 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

add_restraints(){ 
    touch empty.ndx 
    echo "r $RESNAME & a $RESLASTHEAVY
    r GLY & ri 6 & a CA
    r GLY & ri 17 & a CA
    q" > selection.dat 
    
    cat selection.dat | gmx make_ndx -f system.gro \
        -n empty.ndx \
        -o index.ndx >> $logFile 2>> $errFile 
    check index.ndx 
    
    echo '0 1' | gmx mindist -f system.gro \
        -n index.ndx \
        -o gly6_dist.out >> $logFile 2>> $errFile 
    check gly6_dist.out 
    
    C10=`cat gly6_dist.out | awk '{print $2}'`
    CA=`cat gly6_dist.out | awk '{print $3}'`

    echo "[ bonds ]" > distance_restraints.itp 
    printf ";%6s%6s%6s%8s%8s%8s%12s\n" ai aj func b0 kb >> distance_restraints.itp
    printf "%6s%6s%6s%8s%8s%8s%12s\n" $C10 $CA 6 $glyDist $glyRest >> distance_restraints.itp

    echo '0 2' | gmx mindist -f system.gro \
        -n index.ndx \
        -o gly17_dist.out >> $logFile 2>> $errFile 
    check gly17_dist.out 

    C10=`cat gly17_dist.out | awk '{print $2}'`
    CA=`cat gly17_dist.out | awk '{print $3}'`
    printf "%6s%6s%6s%8s%8s%8s%12s\n" $C10 $CA 6 $glyDist $glyRest >> distance_restraints.itp
} 

system_steep(){
    printf "\t\tSystem steep.............................." 
    if [ ! -f System_steep/system_steep.gro ] ; then 
        create_dir System_steep
        
        cp Build_system/system.top System_steep/.
        cp Build_system/system.gro System_steep/.
        cp Build_system/*.itp System_steep/. 
        cd System_steep

        gmx grompp -f $MDP/system_steep_easy.mdp \
            -p system.top \
            -c system.gro \
            -o system_steep_easy.tpr >> $logFile 2>> $errFile 
        check system_steep_easy.tpr 

        gmx mdrun -deffnm system_steep_easy -nt 1 >> $logFile 2>> $errFile 
        check system_steep_easy.gro 

        gmx grompp -f $MDP/system_steep.mdp \
            -p system.top \
            -c system_steep_easy.gro \
            -o system_steep.tpr >> $logFile 2>> $errFile 
        check system_steep.tpr 

        gmx mdrun -deffnm system_steep >> $logFile 2>> $errFile 
        check system_steep.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

system_nvt(){
    printf "\t\tSystem NVT relaxation....................." 
    if [ ! -f System_nvt/system_nvt.gro ] ; then 
        create_dir System_nvt
        
        cp System_steep/system.top System_nvt/.
        cp System_steep/system_steep.gro System_nvt/.
        cp System_steep/*.itp System_nvt/. 
        cd System_nvt

        gmx grompp -f $MDP/system_nvt_relax.mdp \
            -p system.top \
            -c system_steep.gro \
            -o system_nvt.tpr >> $logFile 2>> $errFile 
        check system_nvt.tpr 

        gmx mdrun -deffnm system_nvt >> $logFile 2>> $errFile 
        check system_nvt.gro 
    
        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

system_npt(){
    printf "\t\tSystem NPT relaxation....................." 
    if [ ! -f System_npt/system_npt.gro ] ; then 
        create_dir System_npt
        
        cp System_nvt/system.top System_npt/.
        cp System_nvt/system_nvt.gro System_npt/.
        cp System_nvt/*.itp System_npt/. 
        cd System_npt

        gmx grompp -f $MDP/system_semiisotropic_npt_relax.mdp \
            -p system.top \
            -c system_nvt.gro \
            -o system_npt.tpr >> $logFile 2>> $errFile 
        check system_npt.tpr 

        gmx mdrun -deffnm system_npt >> $logFile 2>> $errFile 
        check system_npt.gro 
    
        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

heating(){
    printf "\t\tHeating solution to 900K.................." 
    if [ ! -f Heating/heating.gro ] ; then 
        create_dir Heating
        
        cp System_npt/system_npt.gro Heating/. 
        cp System_npt/system.top Heating/. 
        cp System_npt/*.itp Heating/. 
        cd Heating

        includeLine=`cat -n system.top | grep '; Include Position restraint file' | awk '{print $1}' | tail -n1 `
        head -n $includeLine system.top > temp.top 
        echo "
#ifdef POSLIGRES 
#include \"posre_$RESNAME.itp\" 
#endif
        " >> temp.top 
        ((includeLine++)) 
        tail -n +$includeLine system.top >> temp.top
        mv temp.top system.top 

        echo "[ position_restraints ]" > posre_$RESNAME.itp 
        echo ";; Pin sulfur atoms to spacing found on SAM" >> posre_$RESNAME.itp 
        for atom in `grep $RESNAME system_npt.gro  | cut -c9- | grep -E "^[ ]*[CNOS]" | cut -c8- | awk '{print $1}'` ; do 
            printf "%6i%6i%10.f%10.f%10.f\n" $atom 1 1000000  1000000 1000000 >> posre_$RESNAME.itp 
        done 

        gmx grompp -f $MDP/prep_heating_sam.mdp \
            -c system_npt.gro \
            -p system.top \
            -o heating.tpr >> $logFile 2>> $errFile  
        check heating.tpr 

        gmx mdrun -deffnm heating >> $logFile 2>> $errFile 
        check heating.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

heated(){
    printf "\t\tNVT simulation at 900K...................." 
    if [ ! -f Heated_nvt/heated_nvt.gro ] ; then 
        create_dir Heated_nvt
        
        cp Heating/heating.gro Heated_nvt/. 
        cp Heating/system.top Heated_nvt/. 
        cp Heating/*.itp Heated_nvt/. 
        cd Heated_nvt

        gmx grompp -f $MDP/prep_heated_sam.mdp \
            -c heating.gro \
            -p system.top \
            -o heated_nvt.tpr >> $logFile 2>> $errFile 
        check heated_nvt.tpr

        gmx mdrun -deffnm heated_nvt >> $logFile 2>> $errFile 
        check heated_nvt.gro 

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

cooling(){
    printf "\t\tCooling solution to 300K.................." 
    if [ ! -f Cooling/cooling.gro ] ; then 
        create_dir Cooling
        
        cp Heated_nvt/heated_nvt.gro Cooling/. 
        cp Heated_nvt/system.top Cooling/. 
        cp Heated_nvt/*.itp Cooling/. 
        cd Cooling

        gmx grompp -f $MDP/prep_cooling_sam.mdp \
            -c heated_nvt.gro \
            -p system.top \
            -o cooling.tpr >> $logFile 2>> $errFile 
        check cooling.tpr

        gmx mdrun -deffnm cooling >> $logFile 2>> $errFile 
        check cooling.gro 

        echo "Protein-H Protein-H" | gmx trjconv -s cooling.tpr \
            -f cooling.gro \
            -o unfolded.pdb \
            -center \
            -ur compact \
            -pbc mol >> $logFile 2>> $errFile

        clean
        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

production(){
    printf "\t\tProduction run............................" 
    if [[ ! -s prep/Cooling/cooling.gro || ! -s prep/System_npt/system_npt.gro ]] ; then 
        echo ; echo 
        echo "ERROR: You must run prep first!" 
        fi 
    check prep/Cooling/cooling.gro prep/System_npt/system_npt.gro 

    if [ ! -d $fold ] ; then mkdir $fold ; fi 
    cd $fold 
    MOLEC=${fold}_${SOL}

    if [ ! -f Production/${MOLEC}_${totSimTime}ns.gro ] ; then 
        printf "\n" 
        create_dir Production
        
        cp ../prep/System_npt/system.top Production/.
        cp ../prep/System_npt/*.itp Production/. 

        if [ $fold = "folded" ] ; then 
            cp ../prep/System_npt/system_npt.gro Production/.
            startStructure=system_npt.gro 
            MOLEC=folded_$SOL
        else 
            cp ../prep/Cooling/cooling.gro Production/.
            startStructure=cooling.gro 
            MOLEC=unfolded_$SOL
            fi 

        cd Production

        if [ ! -f $MOLEC.tpr ] ; then 
            gmx grompp -f $MDP/production_sam.mdp \
                -p system.top \
                -c $startStructure \
                -o $MOLEC.tpr >> $logFile 2>> $errFile 
        fi 
        check $MOLEC.tpr 

        simTime=0
        while [ $simTime -lt $totSimTime ] ; do 
            ((simTime+=50))
            printf "\t\t\t%10i ns....................." $simTime

            if [ ! -f ${MOLEC}_${simTime}ns.gro ] ; then 
                if [ ! -f $simTime.tpr ] ; then 
                    gmx convert-tpr -s $MOLEC.tpr \
                        -until $((simTime*1000)) \
                        -o $simTime.tpr >> $logFile 2>> $errFile 
                    fi 
                check $simTime.tpr 

                if [ -f $MOLEC.cpt ] ; then 
                    ibrun mdrun_mpi -deffnm $MOLEC \
                        -s $simTime.tpr \
                        -cpi $MOLEC.cpt >> $logFile 2>> $errFile  
                else 
                    ibrun mdrun_mpi -deffnm $MOLEC \
                        -s $simTime.tpr >> $logFile 2>> $errFile
                    fi 
                check $MOLEC.gro 
                mv $MOLEC.gro ${MOLEC}_${simTime}ns.gro 
                check ${MOLEC}_${simTime}ns.gro 
                printf "Success\n" 
            else       
                printf "Skipped\n" 
                fi 
            done 

        printf "\n" 
        clean
        cd ../
    else
        printf "Skipped\n"
        fi  
    cd ../
} 

dssp(){
    printf "\t\tRunning dssp analysis....................." 
    if [ ! -f dssp/helen.nrt ] ; then 
        create_dir dssp
        cd dssp
        clean ##clean early. One of the outputs of gmx do_dssp is a *.dat file. We don't want to delete this while cleaning. 

        if [ ! -f scount.xvg ] ; then 
            echo 'Protein' | gmx do_dssp -f ../Production/$MOLEC.xtc \
                -s ../Production/$MOLEC.tpr \
                -ver 1 \
                -sss HGI \
                -ssdump ssdump.dat \
                -o ss.xpm \
                -a area.xpm \
                -ta totarea.xvg \
                -aa averarea.xvg \
                -sc scount.xvg >> $logFile 2>> $errFile
        fi 
        check scount.xvg ss.xpm area.xpm 

        if [ ! -f area.eps ] ; then 
            gmx xpm2ps -f area.xpm \
                -di $TOP/m2p_files/ps.m2p \
                -by 10 \
                -o area.eps >> $logFile 2>> $errFile 
        fi 
        check area.eps

        if [ ! -f area.png ] ; then 
            ps2pdf area.eps 
            convert area.pdf area.png 
        fi 

        if [ ! -f ss.eps ] ; then 
            gmx xpm2ps -f ss.xpm \
                -di $TOP/m2p_files/ps.m2p \
                -by 10 \
                -o ss.eps >> $logFile 2>> $errFile 
        fi 
        check ss.eps

        if [ ! -f ss.png ] ; then 
            ps2pdf ss.eps 
            convert ss.pdf ss.png 
        fi 

        ##Cut out helen.nrt from scount.xvg. 
        echo "#!/usr/bin/env python
import numpy as np  

with open('scount.xvg') as f : 
    filelines = f.readlines() 

header=0
for line in filelines : 
    if line.startswith('#') or line.startswith('@') : 
        header += 1
        if 'A-Helix' in line : 
            alpha=line.split()[1][1:] 
        if '3-Helix' in line : 
            three=line.split()[1][1:] 
        if '5-Helix' in line : 
            five =line.split()[1][1:] 
        if 'B-Sheet' in line : 
            bsheet = line.split()[1][1:] 
        if 'B-Bridge' in line : 
            bbridge= line.split()[1][1:] 
header -= 2

data = np.genfromtxt('scount.xvg',skip_header=header)
col1 = data[:,0]

try : 
    col2a = data[:,int(alpha)+1]
except NameError : 
    col2a = np.zeros(len(data[:,0])) 
try : 
    col2b = data[:,int(three)+1]
except NameError : 
    col2b = np.zeros(len(data[:,0])) 
try : 
    col2c = data[:,int(five)+1]
except NameError : 
    col2c = np.zeros(len(data[:,0])) 
col2=col2a + col2b + col2c

try : 
    col3a = data[:,int(bsheet)+1]
except NameError : 
    col3a = np.zeros(len(data[:,0]))  
try : 
    col3b = data[:,int(bbridge)+1]
except NameError : 
    col3b = np.zeros(len(data[:,0]))  
col3=col3a + col3b

col4 = 18 - col2 - col3

for i in range(len(col1)) : 
    print \"%5i%5i%5i%5i\"%(col1[i],col2[i],col3[i],col4[i])" > cut_helen.py 

        python cut_helen.py > helen.nrt 
        check helen.nrt 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

rgyr(){
    printf "\t\tCalculating radius of gyration............" 
    if [ ! -f rgyr/gyrate.xvg ] ; then 
        create_dir rgyr
        cd rgyr
        clean 

        echo 'Protein' | gmx gyrate -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -o gyrate.xvg >> $logFile 2>> $errFile 
        check gyrate.xvg

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

minimage(){
    printf "\t\tCalculating minimum image................." 
    if [ ! -f minimage/mindist.xvg ] ; then 
        create_dir minimage
        cd minimage
        clean 
        
        echo 'Protein' | gmx mindist -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -pi \
            -od mindist.xvg >> $logFile 2>> $errFile 
        check mindist.xvg 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
} 

rdf(){
    printf "\t\tCalculating RDFs.........................." 
    if [ ! -f rdf/sam_lys.xvg ] ; then 
        create_dir rdf
        cd rdf
        clean 

        gmx rdf -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -ref 'name CD1 CD2 and resname LEU' \
            -sel 'resname SOL and name OW' \
            -o leu_wat.xvg >> $logFile 2>> $errFile 
        check leu_wat.xvg 

        gmx rdf -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -ref 'name NZ and resname LYS' \
            -sel 'resname SOL and name OW' \
            -o lys_wat.xvg >> $logFile 2>> $errFile 
        check lys_wat.xvg 

        gmx rdf -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -ref 'name CD1 CD2 and resname LEU' \
            -sel 'name C10 and resname LIG' \
            -o sam_leu.xvg >> $logFile 2>> $errFile 
        check sam_leu.xvg 

        gmx rdf -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -ref 'name NZ and resname LYS' \
            -sel 'name C10 and resname LIG' \
            -o sam_lys.xvg >> $logFile 2>> $errFile 
        check sam_lys.xvg 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

nopbc(){
    printf "\t\tCreating nopbc copy......................."
    if [ ! -f nopbc/fit.xtc ] ; then
        create_dir nopbc
        cd nopbc
        clean

        if [ $fold = "folded" ] ; then
            startStructure=system_npt.gro
            MOLEC=folded_$SOL
        else
            startStructure=cooling.gro
            MOLEC=unfolded_$SOL
            fi

        echo "Protein Protein" | gmx trjconv -f ../Production/$startStructure \
            -s ../Production/$MOLEC.tpr \
            -pbc mol \
            -ur compact \
            -center \
            -o nopbc.gro >> $logFile 2>> $errFile
        check nopbc.gro

        echo "Protein Protein" | gmx trjconv -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -pbc mol \
            -ur compact \
            -center \
            -o nopbc.xtc >> $logFile 2>> $errFile
        check nopbc.xtc

        echo "1 | 14" > selection.dat
        echo "q" >> selection.dat

        cat selection.dat | gmx make_ndx -f ../Production/$MOLEC.tpr >> $logFile 2>> $errFile
        check index.ndx

        echo "Protein_LIG" | gmx trjconv -f ../Production/$startStructure \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -pbc mol \
            -ur compact \
            -o fit.gro >> $logFile 2>> $errFile
        check fit.gro

        echo "Protein Protein_LIG" | gmx trjconv -f ../Production/$MOLEC.xtc \
            -s ../Production/$MOLEC.tpr \
            -n index.ndx \
            -fit transxy \
            -pbc res \
            -dt 40 \
            -ur compact \
            -o fit.xtc >> $logFile 2>> $errFile
        check fit.xtc

        printf "Success\n"
        cd ../
    else
        printf "Skipped\n"
        fi
}

rmsd(){
    printf "\t\tCalculating RMSD.........................."
    if [ ! -f rmsd/rmsd.xvg ] ; then
        create_dir rmsd
        cd rmsd
        clean

        echo "Backbone Backbone" | gmx rms -f ../nopbc/nopbc.xtc \
            -s $fileName \
            -o rmsd.xvg >> $logFile 2>> $errFile
        check rmsd.xvg

        printf "Success\n"
        cd ../
    else
        printf "Skipped\n"
        fi
}

cd_spectra(){
    printf "\t\tExtracting frames for DichroCalc.........." 
    if [ ! -f cd_spectra/extract.tar.gz ] ; then 
        create_dir cd_spectra
        cd cd_spectra
        clean 
        create_dir extract
        cd extract
        clean 

        echo 'Protein' | gmx trjconv -f ../../Production/$MOLEC.xtc \
            -s ../../Production/$MOLEC.tpr \
            -dt 100 \
            -pbc mol \
            -sep \
            -o frame.pdb >> $logFile 2>> $errFile 
        cd ../

        tar cvfz extract.tar.gz extract >> $logFile 2>> $errFile 
        check extract.tar.gz 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

cluster(){
    printf "\t\tCluster analysis.........................." 
    if [ ! -f cluster/clust-size.xvg ] ; then 
        create_dir cluster
        cd cluster
        clean 

        echo "Backbone Backbone" | gmx rms -s ../Production/$MOLEC.tpr \
            -f ../Production/$MOLEC.xtc \
            -dt 100 \
            -m rmsd.xpm >> $logFile 2>> $errFile 
        check rmsd.xpm 

        echo "Backbone Protein" | gmx cluster -s ../Production/$MOLEC.tpr \
            -f ../Production/$MOLEC.xtc \
            -dm rmsd.xpm \
            -cutoff 0.25 \
            -method gromos \
            -dt 100 \
            -tr clust-trans.xpm \
            -cl clusters.pdb \
            -sz clust-size.xvg >> $logFile 2>> $errFile 
        check clust-size.xvg clusters.pdb 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

good-turing(){
    printf "\t\tGood-Turing Stats........................." 
    gtScripts=$TOP/good-turing_scripts
    if [[ ! -f good-turing/good_turing.rmsd.tar && -f cluster/rmsd.xpm ]] ; then 
        create_dir good-turing
        cd good-turing
        clean 

        python $gtScripts/xpm2dat.py ../cluster/rmsd.xpm 
        check rmsd.dat 

        ##https://github.com/pkoukos/GoodTuringMD
        ##  Adapted to avoid interacitve file name input
        echo "source('$gtScripts/Good_Turing.R')" > run_gt_stat.R

        RScript --no-save --no-restore --verbose run_gt_stat.R > good-turing.log 2>&1 
        check good_turing.rmsd.tar

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

order(){
    printf "\t\tCalculating order of SAM.................." 
    if [ ! -f order/density.xvg ] ; then 
        create_dir order
        cd order
        clean 

        echo "r LIG & a S1 C1 C2 C3 C4 C5 C6 C7 C8 C9 C10" > selection.dat 
        echo "q" >> selection.dat 

        touch empty.ndx 
        cat selection.dat | gmx make_ndx -f ../Production/$MOLEC.tpr \
            -n empty.ndx \
            -o heavyLIG.ndx >> $logFile 2>> $errFile 
        check heavyLIG.ndx 

        gmx density -s ../Production/$MOLEC.tpr \
            -f ../Production/$MOLEC.xtc \
            -n heavyLIG.ndx \
            -d Z \
            -sl 1000 \
            -e 1000 \
            -o density.xvg >> $logFile 2>> $errFile 
        check density.xvg

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

density(){
    equilTime=150  ## ns 
    printf "\t\tCalculating density profile of peptide...." 
    if [[ ! -f denisty/water.xvg || ! -f density/density.xvg || ! -f density/densmap.dat ]] ; then 
        create_dir density
        cd density
        clean 

        #To avoid articfact in the rounding of positions to nearest bin, 
        #   (ie, it seem gmx density only rounds up, so peaks appear when rounding up is more likely)
        #   we need to set the number of bins roughly equal to to the precision of the xtc file. 
        zdim=`tail -n 1 ../Production/system_npt.gro | awk '{print $3}'`
        sliceNum=`echo "$zdim * 1000" | bc -l | awk '{printf "%i", $0}' `
        if [ ! -f water.xvg ] ; then 
            echo 'Water' | gmx density -s ../Production/$MOLEC.tpr \
                -f ../Production/$MOLEC.xtc \
                -b $((equilTime*1000)) \
                -d Z \
                -sl $sliceNum \
                -dens number \
                -o water.xvg >> $logFile 2>> $errFile 
        fi
        check water.xvg

        if [ ! -f density.xvg ] ; then 
            echo 'Protein' | gmx density -s ../Production/$MOLEC.tpr \
                -f ../Production/$MOLEC.xtc \
                -b $((equilTime*1000)) \
                -d Z \
                -sl $sliceNum \
                -dens number \
                -o density.xvg >> $logFile 2>> $errFile 
        fi
        check density.xvg

        if [ ! -f densmap.dat ] ; then 
            if [ ! -f nopbc.xtc ] ; then 
                echo 'Protein Protein' | gmx trjconv -s ../Production/$MOLEC.tpr \
                    -f ../Production/$MOLEC.xtc \
                    -pbc mol \
                    -o nopbc.xtc >> $logFile 2>> $errFile 
            fi 
            check nopbc.xtc 

            if [ ! -f nopbc.gro ] ; then 
                echo 'Protein Protein' | gmx trjconv -s ../Production/$MOLEC.tpr \
                    -f ../Production/system_npt.gro \
                    -pbc mol \
                    -o nopbc.gro >> $logFile 2>> $errFile 
            fi 
            check nopbc.gro 

            if [ ! -f fit.xtc ] ; then 
                echo 'Protein Protein' | gmx trjconv -s nopbc.gro \
                    -f nopbc.xtc \
                    -fit rotxy+transxy \
                    -o fit.xtc >> $logFile 2>> $errFile 
            fi 
            check fit.xtc 

            if [ ! -f densmap.xpm ] ; then 
                echo 'Protein' | gmx densmap -s nopbc.gro \
                    -f fit.xtc \
                    -b $((equilTime*1000)) \
                    -unit count \
                    -o densmap.xpm >> $logFile 2>> $errFile 
            fi 
            check densmap.xpm 

            $TOP/analysis_scripts/xpm2dat.py densmap.xpm 
        fi
        check densmap.dat 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

hbond(){
    equilTime=150  ## ns 
    printf "\t\tAnalyzing hbonds to peptide..............." 
    if [[ ! -f hbond/sidechain.xvg || ! -f hbond/mainchain.xvg || -f hbond/protein.xvg ]] ; then 
        create_dir hbond
        cd hbond
        clean 

        if [ ! -f protein.xvg ] ; then 
            echo 'Protein Water' | gmx hbond -s ../Production/$MOLEC.tpr \
                -f ../Production/$MOLEC.xtc \
                -b $((equilTime*1000)) \
                -num protein.xvg >> $logFile 2>> $errFile 
        fi 
        check protein.xvg 

        if [ ! -f mainchain.xvg ] ; then 
            echo 'MainChain+H Water' | gmx hbond -s ../Production/$MOLEC.tpr \
                -f ../Production/$MOLEC.xtc \
                -b $((equilTime*1000)) \
                -num mainchain.xvg >> $logFile 2>> $errFile 
        fi 
        check mainchain.xvg 

        if [ ! -f sidechain.xvg ] ; then 
            echo 'SideChain Water' | gmx hbond -s ../Production/$MOLEC.tpr \
                -f ../Production/$MOLEC.xtc \
                -b $((equilTime*1000)) \
                -num sidechain.xvg >> $logFile 2>> $errFile 
        fi
        check sidechain.xvg 

        printf "Success\n" 
        cd ../
    else
        printf "Skipped\n"
        fi  
}

main
