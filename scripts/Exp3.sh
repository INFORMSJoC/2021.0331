#Script used to run Experiment 3 of the paper using a linux-based  batch processing system:

#!/bin/bash
#SBATCH -n 8
#SBATCH -t 72:00:00
#SBATCH -o Rou-Exp3.%j.out
#SBATCH -e Rou-Exp3.%j.err

It=30
n=100
S=999
T=-123

for ((k=1; k <= 10; k++))
do	
	for ((i=1; i <= $It; i++))
	do

	./CR_ROU -n $n -s $S $T -r -100 100 -c 1 -f exp3.out

	S=$[$S-$i*3]
	T=$[$S+$i*5]

	done
	n=$[$n+100]