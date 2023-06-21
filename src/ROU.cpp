#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h> //For random numbers
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <chrono>
#include <math.h>
#include <vector>
#include <ctime>
#include <gmp.h>
#include "GF.h"
#include "cmdOpt.h"
#include "Exact_LU_Alg.h"
#include "IPGE.h"
#include "QGE.h"

using namespace std;

struct RunStats{

	//Algorithm running times
	long long ipge_time, ipge_Ah_time, ipge_sub_time, ipge_up_time;
	//To store the ROU special case 2 outputs: row and column indices used to generate v (i.e., when forcing a subvector of v to be linearly dependent on A's columns) and the related number of times the ROU encountered special case 2, as described in Escobedo (2023) 
	int* sc_outputs; 
	//For verifying the correctness of the algorithm outputs
	bool ipge_check;
	//To print out timing statistics
	char header[200];
};

//See extended definitions for the description of each function
void printStats(cmdOpt&, RunStats&, IPGE&);
void check_sols(cmdOpt&, RunStats&, IPGE&);

int main (int argc, char* argv[]) 
{
	//Initialize cmdOpt data structure; it holds the initialization parameters and user-provided run options. Initialize also struct for calculating run times
	cmdOpt run1;	
	RunStats r_stats;

	//Initialize IPGE objects, which will store the updated REF LU factorizations
	//prob1 will calculate and eventually store the updated REF LU factorization associated with the rank one update algorith from Escobedo (2023); to be denoted as ROU for short
	//prob2 will calculate the REF LU factorization of the updated matrix Ah (i.e., factorization from scratch)
	IPGE prob1, prob2;
	
	//Process command line input to load this run's cmdOpt
	if(run1.process_cmd_args(argc, argv)==0)
	{
		run1.show_usage(argv[0]);
		return 0;
	}

	prob1 = IPGE(run1);

	//Save initial matrix A to A0 (needed for ROU) 
	prob1.A0 = GFz_mat_init_zeros(prob1.rows, prob1.cols);
	GFz_mat_copy(prob1.A, prob1.A0, prob1.rows, prob1.cols);		

	//Run REF LU factorization algorithm (and record run time)
	//======================================Timer Start======================================
	auto t1 = std::chrono::high_resolution_clock::now();
	prob1.REF_LU_noPivoting(run1.print);
	auto t2 = std::chrono::high_resolution_clock::now();
	//========================================Timer End======================================
	r_stats.ipge_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();	
	cout << "REF-LU(A) completed...";

	
	//Run the rank-one REF LU update (and record run time); the algorithm returns the column and row indices used to generate v (i.e., to force special case 2 to occur) and the related number of special case 2 calls required by algorithm
	//======================================Timer Start======================================
	t1 = std::chrono::high_resolution_clock::now();
	r_stats.sc_outputs = prob1.Rank_One_Update(run1.print, run1.up_idx_1, run1.up_idx_2);
	t2 = std::chrono::high_resolution_clock::now();
	//========================================Timer End======================================
	r_stats.ipge_up_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	cout << "Rank-one update of REF-LU(A) completed...";
		
	
   //Copy the inputs to probl REF LU factorization in order to calculate Ah and REF-LU(Ah) from scratch
	prob2 = IPGE(run1);
	
	//Copy also the respective row and column index permutations associated with prob1
	for (int i=0; i < prob2.rows; i++)
		prob2.row_idx[i] = prob1.row_idx[i];
	
	for (int j=0; j < prob2.cols; j++)
		prob2.col_idx[j] = prob1.col_idx[j];
	
	//Copy also the update vectors that were used with prob1 
	prob2.v = GFz_vec_init_zeros(prob2.rows);
	prob2.w = GFz_vec_init_zeros(prob2.cols);
	
	GFz_vec_copy(prob1.v, prob2.v, prob2.rows);
	GFz_vec_copy(prob1.w, prob2.w, prob2.cols);

	mpz_t temp;
	mpz_init(temp);
	
	//Calculate the update matrix Ah from its definition (A + vw^T)
	for (int i=0; i < prob2.rows; i++)
	{
		for (int j=0; j < prob2.cols; j++)
		{
			mpz_mul(temp, prob2.v[prob2.row_idx[i]], prob2.w[prob2.col_idx[j]]);
			mpz_add(prob2.A[prob2.row_idx[i]][prob2.col_idx[j]], prob2.A[prob2.row_idx[i]][prob2.col_idx[j]], temp);				
		}			
	}
	
	//Calculate the REF LU factorization of Ah
	//Run REF LU factorization algorithm (and record run time)
	//======================================Timer Start======================================
	t1 = std::chrono::high_resolution_clock::now();
	prob2.REF_LU_noPivoting(run1.print);
	t2 = std::chrono::high_resolution_clock::now();
	//========================================Timer End======================================
	cout << "REF-LU(Ah) completed\n";
	
	r_stats.ipge_Ah_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	
	//Check algorithm correctness, if desired; then print summary run statistics
	if(run1.check_sol)
		check_sols(run1, r_stats, prob1);
	cout << run1;
	printStats(run1, r_stats, prob1);

	return 0;
}

//Prints the results of the experiments to a file; each time the program runs, a new row is appended to the file
void printStats(cmdOpt &run, RunStats &r_stats, IPGE &prob)
{
	ofstream runStats;
	string emptyCell = "\t";

	if (!std::ifstream(run.outFile_name))
		sprintf(r_stats.header, "Dimen\tUp_Idx_1\tUp_Idx_2\tSC2_calls\tSeed1\tSeed2\t[lb,ub]\tREF-LU(A)_s\tROU_s\tREF-LU(Ah)\tSpeedup(A)\tSpeedup(Ah)\tCheck\n");
	else
		sprintf(r_stats.header,	"");

	runStats.open(run.outFile_name, ios::out | ios::app);

	runStats << r_stats.header << run.dim << "\t"
		     << r_stats.sc_outputs[0]+1 << "\t" << r_stats.sc_outputs[1]+1 
			 << "\t" << r_stats.sc_outputs[2] << "\t" 
			 << run.seed1 << "\t" << run.seed2 << "\t" << "[" 
			 << run.lb << "," << run.ub << "]\t";

	runStats << (double)r_stats.ipge_time / 1000 << "\t" 
			 << (double)r_stats.ipge_up_time / 1000 << "\t"
			 << (double)r_stats.ipge_Ah_time / 1000 << "\t"
			 << (double)r_stats.ipge_time/r_stats.ipge_up_time << "\t" 
			 << (double)r_stats.ipge_Ah_time/r_stats.ipge_up_time << "\t" ;

	if (run.check_sol)
		runStats << r_stats.ipge_check << "\t";
	else
		runStats << emptyCell;
	
	runStats << endl;

	cout << "\n\n=========================Run Stats==========================" 
		 << "\nREF_LU_A Time (milliseconds): " << r_stats.ipge_time
		 << "\nREF_Up Time (milliseconds): " << r_stats.ipge_up_time
		 << "\nREF_LU_Ah Time (milliseconds): " << r_stats.ipge_Ah_time
		 << "\n#SC2 Calls Required: " << r_stats.sc_outputs[2]
		 << "\n============================================================\n";
	runStats.close();	
}

//Checks that the output of the ROU algorithm is correct by verifying that the results of REF forward substitution (REF FS) with REF-LU(A) and with REF-LU(Ah) match on both v and w
void check_sols(cmdOpt &run, RunStats &r_stats, IPGE &prob)
{
	int v_check = 0;
	int w_check = 0;
	
	//Check that REF FS with REF-LU(A) on v yields the same result as REF FS with REF-LU(Ah) on v
	cout << "Verification: ";
	
	//Set to 0 by default
	r_stats.ipge_check = 0;
	
	mpz_t *temp_sub1;
	mpz_t *temp_sub2;
	
	temp_sub1 = GFz_vec_init_zeros(prob.rows);
	temp_sub2 = GFz_vec_init_zeros(prob.rows);
	
	temp_sub1 = prob.FSub(prob.A,prob.v,prob.rows,false);
	temp_sub2 = prob.FSub(prob.Ah,prob.v,prob.rows,false);
	
	if(run.print)
	{
		GF_idx_rawPrint(temp_sub1, prob.row_idx);
		GF_idx_rawPrint(temp_sub2, prob.row_idx);
	}
	
	v_check = GFz_vecComp(temp_sub1, temp_sub2, prob.rows);
	
	//Check that REF FS with REF-LU(A^T) on w yields the same result as REF FS with REF-LU(Ah^T) on w
	temp_sub1 = prob.FSub(prob.A,prob.w,prob.rows,true);
	temp_sub2 = prob.FSub(prob.Ah,prob.w,prob.rows,true);
	
	if(run.print)
	{
		GF_idx_rawPrint(temp_sub1, prob.row_idx);
		GF_idx_rawPrint(temp_sub2, prob.row_idx);
	}
	
	w_check = GFz_vecComp(temp_sub1, temp_sub2, prob.rows);
	
	if(v_check == 1 and w_check == 1)
	{
		r_stats.ipge_check = 1;
		cout << "Solutions match";		
	}
	
	else
	{
		r_stats.ipge_check = 0;		
		cout << "Solutions DO NOT match";
	}	
}