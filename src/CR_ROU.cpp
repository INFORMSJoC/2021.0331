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
	long long ipge_time, up_CRU_time, up_ROU_time, ipge_sub_time;
	//To reflect the number of times the ROU encountered special case 2, as described in Escobedo (2023)
	int num_perms; 
	//For verifying the correctness of the algorithm outputs
	bool ipge_check;
	//To print out timing statistics
	char header[200];
};

//See extended definitions for the description of each function
void printStats(cmdOpt&, RunStats&, IPGE&);

void check_CRU_ROU_sols(cmdOpt&, RunStats&, IPGE&, IPGE&);

int main (int argc, char* argv[]) 
{
	//Initialize cmdOpt data structure; it holds the initialization parameters and user-provided run options. Initialize also struct for calculating run times
	cmdOpt run1;	
	RunStats r_stats;

	//Initialize IPGE objects, which will store the updated REF LU factorizations
	//prob1 will calculate and eventually store the updated REF LU factorization associated with the rank one update algorith from Escobedo (2023); to be denoted as ROU for short
	//prob2 will calculate and eventually store the updated REF LU factorization associated with the push-and-swap column replacement algorithm from Escobedo and Moreno-Centeno (2017); to be denoted as CR for short
	IPGE prob1, prob2;

	//Process command line input to load this run's cmdOpt
	if(run1.process_cmd_args(argc, argv)==0)
	{
		run1.show_usage(argv[0]);
		return 0;
	}
	
	prob1 = IPGE(run1);

	//Copy initial matrix A to A0 (needed for ROU) 
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
	
	prob1.ex_idx = 0; //Setting the exiting index (i.e., the replaced column) as the first index, per the specifications of experiment 3 of Escobedo (2023)
		
	//Copy the output REF LU factorization for later use with the CR algorithm (the standard approach for performing a column replacement update)
	prob2 = IPGE(run1);	  
	 
	prob2.A0 = GFz_mat_init_zeros(prob1.rows, prob1.cols);
	GFz_mat_copy(prob1.A, prob2.A, prob1.rows, prob1.cols);	
	GFz_mat_copy(prob1.A0, prob2.A0, prob1.rows, prob1.cols);	
			
	
	//Copy also the respective row and column index permutations associated with prob1
	prob2.row_idx.clear();
	prob2.col_idx.clear();
	
	for (int i = 0; i < prob2.rows; i++)
		prob2.row_idx.push_back(prob1.row_idx[i]);
	
	for (int j = 0; j < prob2.cols; j++)
		prob2.col_idx.push_back(prob1.col_idx[j]);
	
	//REF LU column replacement update as performed via push-and-swap approach 
	//Perform the update (and record run time)
	//======================================Timer Start======================================
	t1 = std::chrono::high_resolution_clock::now();
	prob1.Update_slim();
	t2 = std::chrono::high_resolution_clock::now();
	//========================================Timer End======================================
	r_stats.up_CRU_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	cout << "Column replacement of REF-LU(A) via P&S completed...";
	
	
	//Begin REF LU column replacement update as performed via a REF rank-one update 
	
	vector<int> temp_col_idx;
	for (int j = 0; j < prob1.num_ent; j++)
		temp_col_idx.push_back(j);	
	
	//Calculate the v-vector(s) using A0 and the entering columns (A_plus0) from prob1 
	prob2.num_ent = prob1.num_ent;
	
	prob2.A_plus0 =  GFz_mat_init_zeros(prob2.rows, prob2.num_ent);
	GFz_mat_copy(prob1.A_plus0, prob2.A_plus0, prob2.rows, prob2.num_ent);	
	
	for (int i=0; i < prob2.rows; i++)					
		for (int j = 0; j < prob2.num_ent; j++)
			mpz_sub (prob2.A_plus0[prob2.row_idx[i]][j], prob2.A_plus0[prob2.row_idx[i]][j], prob2.A0[prob2.row_idx[i]][j]);					
		
	
	//Run the algorithm (and record run time)
	//======================================Timer Start======================================
	t1 = std::chrono::high_resolution_clock::now();
	r_stats.num_perms = prob2.Rank_One_Update_CR(run1.print,0);
	t2 = std::chrono::high_resolution_clock::now();
	//========================================Timer End======================================
	r_stats.up_ROU_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	cout << "Column replacement of REF-LU(A) via ROU completed\n";

	//Check algorithm correctness, if desired; then print summary run statistics
	if(run1.check_sol)
		check_CRU_ROU_sols(run1, r_stats, prob1, prob2);
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
		sprintf(r_stats.header, "Dimen\tSC2_calls\tSeed1\tSeed2\t[lb,ub]\tREF-LU(A)_s\tCRU_s\tROU_s\tSpeedup(REF-LU/CRU)\tSpeedup(REF-LU/ROU)\tSpeedup(ROU/CRU)\tCheck\n");
	else
		sprintf(r_stats.header,	"");

	runStats.open(run.outFile_name, ios::out | ios::app);

	runStats << r_stats.header << run.dim 
			 << "\t" << r_stats.num_perms << "\t" 
			 << run.seed1 << "\t" << run.seed2 << "\t" << "[" 
			 << run.lb << "," << run.ub << "]\t";

	runStats << (double)r_stats.ipge_time / 1000 << "\t" 
			 << (double)r_stats.up_CRU_time / 1000 << "\t"
			 << (double)r_stats.up_ROU_time / 1000 << "\t"
			 << (double)r_stats.ipge_time/r_stats.up_CRU_time << "\t" 
			 << (double)r_stats.ipge_time/r_stats.up_ROU_time << "\t"
			 << (double)r_stats.up_ROU_time/r_stats.up_CRU_time << "\t"; 

	if (run.check_sol)
		runStats << r_stats.ipge_check << "\t";
	else
		runStats << emptyCell;
	
	runStats << endl;

	cout << "\n\n=========================Run Stats==========================" 
		 << "\nREF_LU_A Time (milliseconds): " << r_stats.ipge_time
		 << "\nREF_CRU Time (milliseconds): " << r_stats.up_CRU_time
		 << "\nREF_ROU Time (milliseconds): " << r_stats.up_ROU_time
		 << "\n#SC2 Calls Required: " << r_stats.num_perms
		 << "\n============================================================\n";

	runStats.close();	
}

//Checks the output of both tested approaches for performing a column replacement update of REF LU
void check_CRU_ROU_sols(cmdOpt &run, RunStats &r_stats, IPGE &prob1, IPGE &prob2)
{
	//Check that REF forward substitution with REF-LU(A) on v yields the same result as REF forward substitution with REF-LU(Ah) on v
	cout << "Verification: ";
	
	//Set to 0 by default
	r_stats.ipge_check = 0;
	
	prob1.set_RHS(run.density);
			
	//======================================Timer Start======================================
	auto t1 = std::chrono::high_resolution_clock::now();
	prob1.FSub();			
	prob1.BSub();
	auto t2 = std::chrono::high_resolution_clock::now();
	//========================================Timer End======================================
	r_stats.ipge_sub_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	
	prob1.set_IPGE_Bq_Xq();
	
	GFz_mat_copy(prob2.Ah, prob2.A, prob2.rows, prob2.cols);	
	
	prob2.set_RHS(run.density);
	
	//======================================Timer Start======================================
	t1 = std::chrono::high_resolution_clock::now();
	prob2.FSub();			
	prob2.BSub();
	t2 = std::chrono::high_resolution_clock::now();
	//========================================Timer End======================================
	r_stats.ipge_sub_time = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
	
	prob2.set_IPGE_Bq_Xq();
	
	//Compare the difference between the two forward substitution vectors; if they are not identical, ipge_check gets set to 0
	mpq_t total;
	mpq_init(total);
	
	mpq_t temp;
	mpq_init(temp);
	
	for(int j=0; j < prob1.cols_B; j++)
	{
		for (int i=0; i < prob1.rows; i++)
		{					
			mpq_sub(temp, prob1.Xq[prob1.row_idx[i]][prob1.col_idx_B[j]], prob2.Xq[prob1.col_idx[i]][prob1.col_idx_B[j]]);	
			mpq_abs(temp, temp);							
			mpq_add(total, total, temp);
		} 
	}
	
	if(mpq_sgn(total)==0)
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