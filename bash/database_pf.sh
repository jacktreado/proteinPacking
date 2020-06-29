#!/bin/bash

# input
db=$1
NMCPTS=$2
partition=$3
time=$4
max=$5

# directories
pdbDir=~/project/pdb/$db
datDir=$pdbDir/dat
gitDir=~/pdb/proteinPacking
pfDir=$pdbDir/pf

# make directories to use
mkdir -p bin
mkdir -p tasks
mkdir -p slurm
mkdir -p out
mkdir -p $pfDir

# compile into binary
binf=bin/"$db"_pf.o
srcloc=$gitDir/src
voroloc=$srcloc/voro++/src
srcf=$gitDir/main/pf_calc_main.cpp

# compile
rm -f $binf
g++ --std=c++11 -I $srcloc -I $voroloc $srcf $srcloc/*.cpp $voroloc/voro++.cc -o $binf 

if [[ ! -f $binf ]]
then
    echo -- binary file does not exist, compile failed.
    exit 1
fi

# create task file
taskf=tasks/"$db"_pf.task
rm -f $taskf

# loop over files
let fcount=0
let arraynum=0
dir=$datDir/*_H.dat


# LOOP OVER FILES. NOTE: will get all files with seed >= startseed but <= end seed
for f in $dir; do
    # increment number of files                                                                                                                                                                                                            
    let fcount=$fcount+1

    # test number of files found, to break up job submissions
    if [[ $fcount -gt $max ]]
    then
	   break
    else
	   let arraynum=$arraynum+1
    fi

    # parse file name
    file=${f##*/}
    baseid=${file%%.dat}

    # get input file
    inputFileStr=$datDir/$file

    # create output file
    outputFileStr=$pfDir/"$baseid"_pf.dat

    # echo command to task file
    echo running params: $inputFileStr $outputFileStr $NMCPTS
    echo cd `pwd` \; ./$binf $inputFileStr $outputFileStr $NMCPTS >> $taskf
done

if [[ ! -f "$taskf" ]]
then
    echo task file not created, ending before job submission
    exit 1
fi

echo -- total number of array runs = $arraynum

# setup slurm files
slurmf=slurm/"$db"_pf.slurm
job_name="$db"_pf
runout=out/"$db"_pf-%a.out
rm -f $slurmf

# echo about time
echo -- running time = $time for $partition

echo -- PRINTING SLURM FILE...
echo \#\!/bin/bash >> $slurmf
echo \#SBATCH --cpus-per-task=1 >> $slurmf
echo \#SBATCH --array=1-$arraynum >> $slurmf
echo \#SBATCH -n 1 >> $slurmf
echo \#SBATCH -p $partition >> $slurmf
echo \#SBATCH -J $job_name >> $slurmf
echo \#SBATCH -o $runout >> $slurmf
echo sed -n \"\$\{SLURM_ARRAY_TASK_ID\}p\" "$taskf" \| /bin/bash >> $slurmf
cat $slurmf

# run sbatch file
echo -- running on slurm in partition $partition
sbatch -t $time $slurmf


# ====================
#       INPUTS
# ====================
# 1. data base (use syntax from directory name)
# 2. NMCPTS (number of monte carlo pts for mass calc)
# 3. partition
# 4. time
# 5. max number of files considered
