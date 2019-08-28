#!/bin/bash

all=false
prep=false
folded=false
unfolded=false
stampede=false
simTime=50

HELP(){
    echo 
    echo "This programs is a wrapper for the run_sam.sh, run_water.sh, and run_tert.sh scripts." 
    echo 
    echo "Usage: $0 [options] " 
    echo "    -s Run on stampede using sbatch " 
#    echo "    -a All systems " 
    echo "    -p All prep    " 
    echo "    -f All three folded " 
    echo "    -u All three unfolded system " 
    echo "    -t Simulation time (Default = 50ns) " 
    echo "    -h Print this help message and exit " 
    echo 
    exit
}

while getopts sapfut:h opt ; do 
    case $opt in 
       s) 
         stampede=true 
         ;; 
#       a) 
#         all=true
#         ;; 
       p) 
         prep=true
         ;; 
       f)
         folded=true
         ;; 
       u) 
         unfolded=true
         ;; 
       t) 
         simTime=$OPTARG
         ;; 
       :) 
         echo " option $OPTARG requires an argument"
         usage
         ;; 
       h) 
         HELP 
         ;; 
       \?) 
         echo "Inavalid option -$OPTARG " 
         HELP 
         ;; 
     esac 
     done 

if ! $prep && ! $folded && ! $unfolded ; then 
    echo "Nothing to do! Exitting. " 
    echo "   Use $0 -h for more information" 
    exit 
    fi 

if [ ! -z $WORK ] && ! $stampede ; then 
    echo "Must use -s flag to run on stampede"
    exit
    fi 
if [ $simTime -lt 50 ] ; then 
    echo "Simulation time must be greater than 50!" 
    exit 
    fi 

solvent="sam"
sol="sam"
samList="LIG C11"
molList="folded unfolded" 

if $prep ; then 
    for sam in $samList ; do 

        printf "\n\n\t\t*** $sam ***\n\n" 

        if ! $stampede ; then 
            bash run_${solvent}.sh -f prep -L $sam
        else 
            echo "#!/bin/bash" > submit_prep_${sam}
            echo >> submit_prep_${sam}
            echo "#SBATCH -J prep_${sam} " >> submit_prep_${sam}
            echo "#SBATCH -o prep_${sam}.o%j" >> submit_prep_${sam}
            echo "#SBATCH -N 1" >> submit_prep_${sam}
            echo "#SBATCH -n 48 " >> submit_prep_${sam}
            echo "#SBATCH -p skx-normal " >> submit_prep_${sam}
            echo "#SBATCH -t 05:00:00" >> submit_prep_${sam}
            echo "#SBATCH -A Ras" >> submit_prep_${sam}
            echo "#SBATCH --mail-user=jeremy_first@utexas.edu" >> submit_prep_${sam}
            echo "#SBATCH --mail-type=all" >> submit_prep_${sam}
            
            echo >> submit_prep_${sam}
            echo "module load gromacs " >> submit_prep_${sam} 
                                                                           
            echo >> submit_prep_${sam}
            echo "bash run_${solvent}.sh -f prep -L $sam" >> submit_prep_${sam}

            sbatch submit_prep_${sam} 
            fi 
        done 
    fi 

if $folded ; then 
    for sam in $samList ; do 

        printf "\n\n\t\t*** $sam ***\n\n" 

        if ! $stampede ; then 
            bash run_${solvent}.sh -f folded -L $sam
        else 
            echo "#!/bin/bash" > submit_fol_${sam}
            echo >> submit_fol_${sam}
            echo "#SBATCH -J fol_${sam} " >> submit_fol_${sam}
            echo "#SBATCH -o fol_${sam}.o%j" >> submit_fol_${sam}
            echo "#SBATCH -N 1" >> submit_fol_${sam}
            echo "#SBATCH -n 48 " >> submit_fol_${sam}
            echo "#SBATCH -p skx-normal " >> submit_fol_${sam}
            echo "#SBATCH -t 48:00:00" >> submit_fol_${sam}
            echo "#SBATCH -A Ras" >> submit_fol_${sam}
            echo "#SBATCH --mail-user=jeremy_first@utexas.edu" >> submit_fol_${sam}
            echo "#SBATCH --mail-type=all" >> submit_fol_${sam}
            
            echo >> submit_fol_${sam}
            echo "module load gromacs " >> submit_fol_${sam} 
                                                                           
            echo >> submit_fol_${sam}
            echo "bash run_${solvent}.sh -f folded -L $sam -t $simTime " >> submit_fol_${sam}

            sbatch submit_fol_${sam} 
            fi 
        done 
    fi 

if $unfolded ; then 
    for sam in $samList ; do 

        printf "\n\n\t\t*** $sam ***\n\n" 

        if ! $stampede ; then 
            bash run_${solvent}.sh -f unfolded -L $sam
        else 
            echo "#!/bin/bash" > submit_unf_${sam}
            echo >> submit_unf_${sam}
            echo "#SBATCH -J unf_${sam} " >> submit_unf_${sam}
            echo "#SBATCH -o unf_${sam}.o%j" >> submit_unf_${sam}
            echo "#SBATCH -N 1" >> submit_unf_${sam}
            echo "#SBATCH -n 48 " >> submit_unf_${sam}
            echo "#SBATCH -p skx-normal " >> submit_unf_${sam}
            echo "#SBATCH -t 48:00:00" >> submit_unf_${sam}
            echo "#SBATCH -A Ras" >> submit_unf_${sam}
            echo "#SBATCH --mail-user=jeremy_first@utexas.edu" >> submit_unf_${sam}
            echo "#SBATCH --mail-type=all" >> submit_unf_${sam}
            
            echo >> submit_unf_${sam}
            echo "module load gromacs " >> submit_unf_${sam} 
                                                                           
            echo >> submit_unf_${sam}
            echo "bash run_${solvent}.sh -f unfolded -L $sam -t $simTime" >> submit_unf_${sam}

            sbatch submit_unf_${sam} 
            fi 
        done 
    fi 

