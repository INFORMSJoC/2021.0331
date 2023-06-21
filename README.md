This code is used to experiment with the REF LU Factorization's rank-one update algorithms introduced in Escobedo (2023); if using this software, please cite this code, whose separate DOI is to be provided in the published journal article. 
In addition, the associated paper should be cited as:

@article{escobedo2023exact,
  title={Exact Matrix Factorization Updates for Nonlinear Programming},
  author={Escobedo, Adolfo R},
  journal={to appear in INFORMS Journal on Computing},
  url={https://arxiv.org/abs/2202.00520},
  year={2023}
}

*****************************************************************Typing <make> in this folder creates two executables*****************************************************************

1) ROU.exe: This executable is used to (a) build the REF LU factorization of a dense matrix A, (b) perform a rank-one update of REF-LU(A), and (3) build the REF LU factorization of the updated matrix, Ah, from scratch for comparison

2) CR_ROU.exe: This executable is used to (a) build the REF LU factorization of a dense matrix A, (b) perform a column replacement update of REF-LU(A) using the REF rank-one update algorithm, and (3) perform a column replacement update of REF-LU(A) using the push-and-swap approach of Escobedo and Moreno-Centeno (2017)

A description of each executable is given below:

1) ROU.exe:
usage: ./ROU

followed by the following options:
    
	-h,--help
        -c,--check_sol  [1,0]           >default=0<  						(verify that the algorithm outputs are correct - requires additional computations)
        -f,--file_out   [string]        >default_output.out<					(specify name of output file where the results of each run should be recorded - recommend >Runs_ROU.out<)
        -i,--upd_idx    [int] [int]     >default=-1,-1< 					(specify indices used to generate v, when forcing part of it to be linearly dependent on the input matrix's columns)
                        									(if upd_idx_1 =-1, the update vector v is not forced to be linearly dependent)
                        									(if upd_idx_1 = 0, v is forced to be linearly dependent using random indices)
                        									(if upd_idx_1 > 0, v is forced to be linearly dependent using first up_idx1 columns and first up_idx_2 rows of original matrix)

        -n,--num_vars   [int]           >default=5<						(Dimension of the square input matrix)
        -p,--print      [1,0]           >default=0<						(Prints to screen the working matrix after each iteration; for best results, use only with n <= 10)
        -r,--range      [int] [int]     >default=-9,9<						(The lower and upper bound for each generated matrix entry)
        -s,--seeds      [int] [int]     >default=11,1400<					(1st seed for nonzero-entries placement, second for entry values)



2) CR_ROU.exe
usage: ./CR_ROU


	-h,--help
        -c,--check_sol  [1,0]           >default=0<  						(verify that the algorithm outputs are correct - requires additional computations)
        -f,--file_out   [string]        >default_output.out<					(specify name of output file where the results of each run should be recorded - recommend >Runs_CR_ROU.out<)
        -n,--num_vars   [int]           >default=5<						(Dimension of the square input matrix)
        -p,--print      [1,0]           >default=0<						(Prints to screen the working matrix after each iteration; for best results, use only with n <= 10)
        -r,--range      [int] [int]     >default=-9,9<						(The lower and upper bound for each generated matrix entry)
        -s,--seeds      [int] [int]     >default=11,1400<					(1st seed for nonzero-entries placement, second for entry values)



*****************************************************************Replicating Experiment 1 of Escobedo (2023)*****************************************************************

Run executable ROU described above multiple times by setting the respective flags as follows:
  	-n  : with <{16,32,64,128,256,512,1024}> specify each of these integers for 30 repetitions
        -s  : with <int int> (specify two different integer random seeds for each repetitions; this can be done with a script.
	-f  : with <exp_v1.out> (so that all runs are recorded in the same file)
	-r  : with <-100 100> (provides the range for each randomly generated integer matrix entry)
	-c  : with <1> (to certify that the algorithm outputs are correct) 
	 

*****************************************************************Replicating Experiment 2 of Escobedo (2023)*****************************************************************

Run executable ROU described above multiple times by setting the respective flags as follows:
  	-n  : with <{16,32,64,128,256,512,1024}> specify each of these integers for 30 repetitions
        -s  : with <int int> (specify two different integer random seeds for each repetitions; this can be done with a script.
	-i  : with <0 0>, separated by a space (this generates the update vector v to be based partially on a linear combination of some of the input matrix columns)
	-f  : with <exp_v2.out> (so that all runs are recorded in the same file)
	-r  : with <-100 100> (provides the range for each randomly generated integer matrix entry)
	-c  : with <1> (to certify that the algorithm outputs are correct) 

*****************************************************************Replicating Experiment 3 of Escobedo (2023)*****************************************************************

Run executable CR_ROU described above multiple times by setting the respective flags as follows:
  	-n  : <{100,200,300,400,500,600,700,800,900,1000}> specify each of these integers for 30 repetitions
        -s  : <int int> (specify two different integer random seeds for each repetitions; this can be done with a script.
	-f  : with <exp_v3.out> (so that all runs are recorded in the same file)
	-r  : with <-100 100> (provides the range for each randomly generated integer matrix entry)
	-c  : with <1> (to certify that the algorithm outputs are correct) 
	 

*****************************************************************
Note that the above flags provide more flexibility in how the algorithm can be tested than what is done in Escobedo (2023).
