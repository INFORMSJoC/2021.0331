#Script used to run Experiment 1 of the paper using a linux-based  batch processing system:

#!/bin/bash
#SBATCH -n 8
#SBATCH -t 72:00:00
#SBATCH -o Rou-Exp1.%j.out
#SBATCH -e Rou-Exp1.%j.err

It=30
n=16
S=-2878
T=-2728

cd ../src/

for ((k=1; k <= 6; k++))
do	
	for ((i=1; i <= $It; i++))
	do

	./ROU -a 0 -n $n -s $S $T -r -100 100 -c 1 -f ../results/exp1.out

	S=$[$S-$i*3]
	T=$[$S+$i*5]

	done
	n=$[$n*2]
