[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Exact Matrix Factorization Updates for Nonlinear Programming

The files included in this subdirectory are the computer outputs associated with the three experiments performed in the paper. 

exp1.out corresponds to Experiment 1. The experiment can be replicated by running scripts/Exp1.sh in a Linux-based batch processing system; alternatively, execute ROU.exe multiple times, setting the respective flags as follows:

  	-n  : with <{16,32,64,128,256,512,1024}> specify each of these integers for 30 repetitions
        -s  : with <int int> (specify two different integer random seeds for each repetitions; this can be done with a script.
	-f  : with <exp_v1.out> (so that all runs are recorded in the same file)
	-r  : with <-100 100> (provides the range for each randomly generated integer matrix entry)
	-c  : with <1> (to certify that the algorithm outputs are correct) 

exp2.out corresponds to Experiment 2. The experiment can be replicated by running scripts/Exp2.sh in a Linux-based batch processing system; alternatively, execute ROU.exe multiple times, setting the respective flags as follows:

  	 -n  : with <{16,32,64,128,256,512,1024}> specify each of these integers for 30 repetitions
        -s  : with <int int> (specify two different integer random seeds for each repetitions; this can be done with a script.
	-i  : with <0 0>, separated by a space (this generates the update vector v to be based partially on a linear combination of some of the input matrix columns)
	-f  : with <exp_v2.out> (so that all runs are recorded in the same file)
	-r  : with <-100 100> (provides the range for each randomly generated integer matrix entry)
	-c  : with <1> (to certify that the algorithm outputs are correct) 

exp3.out corresponds to Experiment 3. The experiment can be replicated by running scripts/Exp3.sh in a Linux-based batch processing system; alternatively, execute CR_ROU.exe multiple times, setting the respective flags as follows:

  	-n  : <{100,200,300,400,500,600,700,800,900,1000}> specify each of these integers for 30 repetitions
        -s  : <int int> (specify two different integer random seeds for each repetitions; this can be done with a script.
	-f  : with <exp_v3.out> (so that all runs are recorded in the same file)
	-r  : with <-100 100> (provides the range for each randomly generated integer matrix entry)
	-c  : with <1> (to certify that the algorithm outputs are correct)