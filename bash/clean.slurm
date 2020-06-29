#!/bin/bash

# BASH SCRIPT to "clean" _H pdb files to make them into more manageable data files using the MATLAB
# code cleanProteinDir.
#
# NOTES:
#     1. Make sure database name is correct, needs to be hardcoded into script
#     2. Make sure matlab file and radii.dat are in same working directory as script
#


#SBATCH --partition=scavenge
#SBATCH --time=0-12:00:00
#SBATCH -J clean
#SBATCH -o clean.out

# load MATLAB
module load MATLAB

# data base name
dbname=dun1_0

# input parameters
inDir=~/project/pdb/$dbname/wH
outDir=~/project/pdb/$dbname/dat
radLoc=`pwd`
frmt=_H.pdb

# make output directory if it doesnt exist yet
mkdir -p $outDir

# write matlab string
MCODE="cleanProteinDir('$inDir','$outDir','$radLoc','$frmt'); quit;"

# submit to cluster
echo Now running $MCODE...
matlab -nodisplay -r "$MCODE"
echo finished running $MCODE

