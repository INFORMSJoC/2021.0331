class cmdOpt{

	public: 
		int alg_choice, dim, lb, ub, num_RHS, seed1, seed2, up_idx_1, up_idx_2;
		double density;
		bool print, check_sol;
		string outFile_name;

		cmdOpt(): 
		outFile_name(""),
		num_RHS(0),
		alg_choice(0),
		dim(5),
		seed1(11),
		seed2(1400),
		lb(-9),
		ub(9),
		up_idx_1(-2),
		up_idx_2(-2),
		density(1.0),
		print(0),
		check_sol(0)
		{}	

		// Process command line arguments 
		int process_cmd_args(int numArgs, char* args[]);
		
		// Displays list of flags that can be specified, along with default values for the respective parameters
		void show_usage(string name);
};

// Process command line arguments
int cmdOpt::process_cmd_args(int numArgs, char* args[])
{
	if(numArgs > 1)
	{
		for (int i = 1; i < numArgs; ++i)
		{
			string arg = args[i]; 
	
			// User help parameter
			if (arg == "-h" || arg == "--help")
				return 0;

			// Choose update algorithm - general ROU update algorithm or special column replacement update (via ROU)
			else if (arg == "-a" || arg == "--algorithm")
			{
				alg_choice = atoi(args[++i]);
				if (alg_choice != 0 && alg_choice != 1)
				{
					cerr << "\n*** Error: alg_choice input restricted to 0 or 1 ***\n";
					return 0;
				}
			}

			// Specifies whether validation of output should be performed (requires more calculations)
			else if (arg == "-c" || arg == "--check_sol")
			{
				++i;
				if (atoi(args[i]) != 0 && atoi(args[i]) != 1)
				{
					cerr << "\n*** Error: Check input restricted to 0 or 1 ***\n";
					return 0;
				}
								
				check_sol = (atoi(args[i]) == 1);
				
			}

			// Percentage of non-zeros in the generated data structures
			else if (arg == "-d" || arg == "--density")
			{
				density = atof(args[++i]);
				if (density < 0 || density > 1)
				{
					cerr << "\n*** Error: Density must be in interval [0,1] ***\n";
					return 0;
				}
			}

			// Name of file where run outputs are stored
			else if (arg == "-f" || arg == "--file_out")
			{
				outFile_name = args[++i];
				
				if (outFile_name.size() < 1)
				{
					cerr << "\n*** Error: No file string input ***\n";
					return 0;
				}
			}

			// Specify whether v should be generated to be partially linearly independent (to force special case 2 adjustments)
			else if (arg == "-i" || arg == "--Update_vector_initiation_indices")
			{
				up_idx_1 = atoi(args[++i])-1;
				up_idx_2 = atoi(args[++i])-1;
				
				if (up_idx_1 < -2 || up_idx_1 > dim-1)
				{
					cerr << "\n*** Error: Column index (upd_idx_1) needed for rank-one update vector initiation must lie between -1 and n***\n";
					return 0;
				}
				
				if (up_idx_1 > -1 && (up_idx_2 < up_idx_1 || up_idx_2 > dim-1))
				{
					cerr << "\n*** Error: When specifying a positive value for up_idx1, the row index (up_idx_2) needed for rank-one update vector initiation must lie 1 and " << dim << "***\n";
					return 0;
				}

				if (alg_choice==1 && (up_idx_1 > -2 || up_idx_2 > -2))
				{
					cerr << "\n*** Error: cannot specify values for v to be partially linearly dependent when applying the column replacement update ***\n";
					return 0;
				}
			}

			// Dimension of basis matrix from which factorization is constructed
			else if (arg == "-n" || arg == "--num_vars")
			{
				dim = atoi(args[++i]);
				if (dim < 1)
				{
					cerr << "\n*** Error: Number of variables must be positive ***\n";
					return 0;
				}
			}

			// Specify whether detailed output should be generated (not recommended for n>10)
			else if (arg == "-p" || arg == "--print")
			{
				++i;
				if (atoi(args[i]) != 0 && atoi(args[i]) != 1)
				{
					cerr << "\n*** Error: Print input restricted to 0 or 1 ***\n";
					return 0;
				}
				print = (atoi(args[i]) == 1);
			}

			// Range from which the integer entries are generated
			else if (arg == "-r" || arg == "--range")
			{
				lb = atoi(args[++i]);
				ub = atoi(args[++i]);
				if(lb > ub)
				{
					cerr << "\n*** Error: Range lower bound must be less than or equal to upper bound ***\n";
					return 0;					
				}
			}

			// Random seeds used to generate the input data
			else if (arg == "-s" || arg == "--seeds")
			{
				seed1 = atoi(args[++i]);
				seed2 = atoi(args[++i]);
			}
			
			else return 0;
		}
		cout << endl;
	}
	return 1;
}

// Displays list of flags that can be specified, along with default values for the respective parameters
void cmdOpt::show_usage(string name)
{
	cerr << "\nUsage: " << name << endl
        << "\nOptions:\n"
        << "\t-h,--help\n"
		<< "\t-a,--alg_choice\t[0,1]\t\t>default=0<\t\n"
		<< "\t\t\t(0 corresponds to the general ROU update algorithm)\n"
		<< "\t\t\t(1 corresponds to the special column replacement update algorithm (via ROU))\n"
		<< "\t-c,--check_sol\t[0,1]\t\t>default=0<\t\n"
		<< "\t-f,--file_out \t[string]\t>default_output.out<\t\n" 
		<< "\t-i,--upd_idx \t[int] [int]\t>default=-1,-1< (relevant only for Up_v4)\n" 
		<< "\t\t\t(if upd_idx_1 = -1, the update vector v is not forced to be linearly dependent)\n"
		<< "\t\t\t(if upd_idx_1 = 0, v is forced to be linearly dependent using random indices)\n"
		<< "\t\t\t(if upd_idx_1 > 0, v is forced to be linearly dependent using first up_idx1 columns and first up_idx_2 rows of original matrix)\n\n"
		<< "\t-n,--num_vars \t[int]\t\t>default=5<\n"
		<< "\t-p,--print \t[0,1]\t\t>default=0<\t\n" 
		<< "\t-r,--range \t[int] [int]\t>default=-9,9<\n"
		<< "\t-s,--seeds \t[int] [int]\t>default=11,1400<\n"
        << endl; 
}

//Overloaded << operator for cmdOpt class; shows key command line arguments provided
std::ostream& operator<<(std::ostream &strm, const cmdOpt &run) 
{
	strm << "\n===================Command Options Input===================="
		 << "\ncheck_sol:\t\t" << run.check_sol
		 << "\nalg_choice:\t\t" << run.alg_choice
		 << "\nout_file:\t\t" << run.outFile_name 
		 << "\nupd_idx:\t\t[" << run.up_idx_1+1 << "," << run.up_idx_2+1 << "]"
		 << "\nnum_vars:\t\t" << run.dim
		 << "\nprint:\t\t\t" << run.print
		 << "\nrange:\t\t\t[" << run.lb << "," << run.ub << "]"
		 << "\nRandom seeds:\t\t" << run.seed1 << "," << run.seed2
		 << "\n============================================================\n";
	
	return strm;
}
