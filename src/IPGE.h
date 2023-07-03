#include <sstream>
class IPGE : public Exact_LU_Alg
{
	public:

		// To keep track of sign changes in the matrix frames
		vector<int> frame_sgn;

		// Default Constructor
		IPGE ();

		// Constructor from user-specified command options
		IPGE (cmdOpt &);

		// Matrix file constructor
		IPGE (string);

		/*-----------------------------------------------------------------------------------------
		                        FUNCTIONS
		-----------------------------------------------------------------------------------------*/
		
		
		// REF LU algorithm with no permutations (i.e., does not account for possibility of ZERO-valued pivots)
		void REF_LU_noPivoting(bool);
		
		//Reduced matrix IPGE pivot (i.e., does not affect previous pivot rows/cols)
		void piv_red();
		
		// This function sets a set of righthand side vectors using random values, according to the specified density of non-zeros
		void set_RHS(double);

		/* 	
			REF forward substitution, as described in Escobedo and Moreno-Centeno (2015)
			The algorithm solves the triangular system (L^{-1})y=b in integer-preserving arithmetic onto the existing right-hand side columns
		*/
		void FSub();

		// Similar as the above function, but the matrix to be forwarded is provided as an argument
		void FSub(mpz_t **, int);

		// Similar as the above function, but the matrix to be forwarded is provided as an argument as well as an index within that matrix (i.e., the substitution is performed only on this column) 
		void FSub(mpz_t **, int, int);
		
		// Similar as the above function, but the function takes REF-LU matrix, a single vector, and number of rows as arguments; also has option for transpose forward substitution; returns forward-substituted vector
		mpz_t *  FSub(mpz_t **, mpz_t *, int, bool);
		
		
		/*
			Single FSub step expanding from the 3-argument F-Sub function defined above
			Performs REF Forward Substitution Stepwise Recusion defined in Escobedo (2023) (with the difference that column step_idx of L and previous pivot are not passed as arguments (they are accessed directly)
			The last argument is used to denote whether the FSub is being performed in the transpose sense (i.e., with U^T rather than with L)
		*/
		void FSub_Step(mpz_t **, int, int, int, bool);
		
		/*	
			Single FSub step expanding from the 3-argument F-Sub function defined above
			Performs REF Forward Substitution Stepwise Recusion defined in Escobedo (2023) (with the difference that column step_idx of L and previous pivot are not passed as arguments (they are accessed directly)
			Last argument is used to denote whether the FSub is being performed in the transpose sense (i.e., with U^T rather than with L)
		*/
		mpz_t * FSub_Step_slim(mpz_t *, int, bool);
		
		/*
			REF backward substitution, as described in Escobedo and Moreno-Centeno (2015)
			The algorithm solves the scaled triangular system Udet(A)x=det(A)y in exact rational arithmetic 
		*/
		void BSub();

		// Initializes the right-hand side vectors and solution vectors
		void set_IPGE_Bq_Xq();



		/*-----------------------------------------------------------------------------------------
		    FUNCTIONS associated with the push-and-swap REF column replacement update
		-----------------------------------------------------------------------------------------*/

		// Resets all values of the binary frame_sgn array to 1 
		void reset_frameSigns();
		
		// REF operation for rewinding all entries within a particular frame (and possibly shifting left)
		void backtrack(int, bool, bool);

		// REF operation for performing a change of pivot operation on a particular frame
		void RwSOP(int, int);
		
		// Backtrack subroutine specific to column-push operation
		void CPush_backtrack(int);

		// Backtrack subroutine specific to column-push operation
		void CPush_backtrack_slim(int);

		// RwSOP subroutine specific to column-push operation
		void CPush_RwSOP(int, int);

		// Reduced RwSOP subroutine specific to column-push operation
		void CPush_RwSOP_slim(int, int);

		// Column permutation subroutine specific to column-push operation
		void CPush_CPerm(int);

		// Reduced column permutation subroutine specific to column-push operation
		void CPush_CPerm_slim(int);
		
		// Changes signs of ALL entries in frames k,...,n (includes extension to entering column)
		void CPush_flipSigns(int, int);
		
		// Reduced function for changing signs of ALL entries in frames k,...,n (includes extension to entering column)
		void CPush_flipSigns_slim(int);
		
		// Column-push algorithm 
		void CPush(int);

		// Reduced column-push algorithm 
		void CPush_slim(int);

		// Reduced column-push algorithm 
		void CPush_slim(int, bool);

		// Push-and-swap REF column replacement update
		void Update();

		// Reduced Push-and-swap REF column replacement update
		void Update_slim();		
		
		
		/*-----------------------------------------------------------------------------------------
		    FUNCTIONS Associated with the REF rank-one update 
		-----------------------------------------------------------------------------------------*/	
		
		// Changes signs of ALL entries in frames k,...,n, for rank-one update
		void ROU_flipSigns(int);
		
		// Changes signs of rank-one update vector starting at the specified index
		void ROU_flipSigns_vec(mpz_t *, int);
		
		// RwSOP subroutine specific to rank-one update
		void ROU_RwSOP(int);
		
		// Backtrack subroutine specific to rank-one update
		void ROU_backtrack(int);
		
		// Sets v vector used in ROU; if given FALSE, v is generated randomly
		// If given TRUE, v is generated as a linear combination of columns of A (using either random or user-provided indices)
		mpz_t* ROU_set_update_vec_v(int,int, bool);
		
		// Sets w vector used in ROU; all entries are generate randomly
		mpz_t* ROU_set_update_vec_w();
				
		// Checks that the next divisor needed to obtain entries of the updated REF-LU factorization is nonzero
		bool ROU_checkNextDivisorGood(mpz_t **, mpz_t*, int, bool);
		
		// Recalculates an entry from REF-LU(Ah0) up to the specified iteration k
		void ROU_recomputeAhDiagEntry(mpz_t **, int, int, int);
		
		// Handles the case when the divisor needed to obtain entries of L is zero; it is a subroutine called, as needed, within the ROU algorithm
		void ROU_Special_Case(int, bool);

		/*
			This is the main algorithm discussed in the paper, with adjustments for special case 2 (SC2) included
			This function returns an int array consisting of (1) the column index used for generating v, (2) the row index used for generating v, and (3) the number of SC2 calls required
		*/
		int* Rank_One_Update(bool, int, int);
		
		/*
			This is the special variant of the ROU algorithm for performing a column replacement update discussed in Escobedo (2023), with adjustments for special case 2 (SC2) included 
			The number of SC2 calls required by the algorithm is returned by this function
		*/
		int* Rank_One_Update_CR(bool, int);
};

// Default Constructor
IPGE::IPGE () {
	Exact_LU_Alg();
}

// Constructor from user-specified command options
IPGE::IPGE (cmdOpt &run) 
{
	load_from_cmdOpt(run);
	A = GFz_mat_init_density(rows, cols, get_lb(), get_ub(), get_seed1(), get_seed2(), get_density(), true);
}

// Matrix file constructor
IPGE::IPGE (string f_name)
{
	string line;
		
	ifstream in_file;
	in_file.open(f_name);
		
	long long temp;

	if (in_file.is_open())
	{	
		getline(in_file,line);
		
		istringstream line_ss(line);
		line_ss >> rows;
		line_ss >> cols;

		A = GFz_mat_init_zeros(rows, cols);
			
		for (int i = 0; i < rows; i++)
		{
			getline(in_file,line);
			istringstream line_ss(line);
			for (int j = 0; j < cols; j++)
			{
				line_ss >> temp;
				mpz_set_si(A[i][j], temp);
			}
		}
	}

	in_file.close();

	stepNum = 0;
	pSpace = ceil(log10(get_ub()+1))+1;

	for (int i = 0; i < rows; i++)
	{
		row_idx.push_back(i);
		row_idx_std.push_back(i);
	}

	for (int j = 0; j < cols; j++)
	{
		col_idx.push_back(j);
		col_idx_std.push_back(j);
	}
}

// REF LU algorithm with no permutations (i.e., does not account for possibility of ZERO-valued pivots)
void IPGE::REF_LU_noPivoting(bool printIt)
{
	if(printIt)
		GFz_idx_formPrint(A, row_idx, col_idx, pSpace, 0);
		
	for(int i=1; i < rows; i++)
	{
		piv_red();
		if(printIt)	
			GFz_idx_formPrint(A, row_idx, col_idx, pSpace, i);
	}
}

//Reduced matrix IPGE pivot (i.e., does not affect previous pivot rows/cols)
void IPGE::piv_red()
{
	mpz_t prod1, prod2;
	mpz_inits(prod1, prod2, NULL); 

	for(int i=stepNum+1; i < rows; i++)
	{
		for(int j=stepNum+1; j < cols; j++)
		{
			mpz_mul (prod1, A[row_idx[stepNum]][col_idx[stepNum]], 
								A[row_idx[i]][col_idx[j]]);

			mpz_mul (prod2, A[row_idx[stepNum]][col_idx[j]], 
							A[row_idx[i]][col_idx[stepNum]]);

			mpz_sub (A[row_idx[i]][col_idx[j]], prod1, prod2);

			if(stepNum > 0)			
				mpz_divexact (A[row_idx[i]][col_idx[j]], A[row_idx[i]][col_idx[j]], 
								A[row_idx[stepNum-1]][col_idx[stepNum-1]]);
		}
	}
	stepNum++;
}

// This function sets a set of righthand side vectors using random values, according to the specified density of non-zeros
void IPGE::set_RHS(double density_)
{
	//Initialize elements in B to 0
	B = GFz_mat_init_density(rows, cols_B, get_lb(), get_ub(), get_seed1()*3, get_seed2()*2, density_, false);	
}

/* 	
	REF forward substitution, as described in Escobedo and Moreno-Centeno (2015)
	The algorithm solves the triangular system (L^{-1})y=b in integer-preserving arithmetic onto the existing right-hand side columns
*/

void IPGE::FSub()
{
	mpz_t prod1, prod2;
	mpz_inits(prod1, prod2, NULL);
			
	//Initalize y-vector(s) to all 0s
	Y = GFz_mat_init_zeros(rows, cols_B);
			
	for (int k = 0; k < cols_B; k++)
	{
		mpz_set(Y[row_idx[0]][k], B[row_idx[0]][k]);

		for (int i = 1; i < rows; i++)
		{
			mpz_set(Y[row_idx[i]][k], B[row_idx[i]][k]);
										
			for(int j = 0; j < i; j++)
			{						
				mpz_mul(prod1, A[row_idx[j]][col_idx[j]], Y[row_idx[i]][k]);
				mpz_mul(prod2, A[row_idx[i]][col_idx[j]], Y[row_idx[j]][k]);
				mpz_sub(Y[row_idx[i]][k], prod1, prod2);

				if (j > 0)
					mpz_divexact(Y[row_idx[i]][k], Y[row_idx[i]][k], A[row_idx[j-1]][col_idx[j-1]]);		
			}	
		}
	}
}

// Similar as the above function, but the matrix to be forwarded is provided as an argument 
void IPGE::FSub(mpz_t **mpz_mat, int mat_cols)
{
	mpz_t prod1, prod2;
	mpz_inits(prod1, prod2, NULL);

	for (int k = 0; k < mat_cols; k++)
	{
		for (int i = 1; i < rows; i++)
		{
			for (int j = 0; j < i; j++)
			{
				mpz_mul(prod1, A[row_idx[j]][col_idx[j]], mpz_mat[row_idx[i]][k]);
				mpz_mul(prod2, A[row_idx[i]][col_idx[j]], mpz_mat[row_idx[j]][k]);
				mpz_sub(mpz_mat[row_idx[i]][k], prod1, prod2);

				if (j > 0)
					mpz_divexact(mpz_mat[row_idx[i]][k], mpz_mat[row_idx[i]][k], A[row_idx[j - 1]][col_idx[j - 1]]);
			}
		}
	}
}

// Similar as the above function, but the matrix to be forwarded is provided as an argument as well as an index within that matrix (i.e., the substitution is performed only on this column) 
void IPGE::FSub(mpz_t **mpz_mat, int mat_cols, int idx_RHS_col)
{
	mpz_t prod1, prod2;
	mpz_inits(prod1, prod2, NULL);
	for (int i = 1; i < rows; i++)
	{
		for (int j = 0; j < i; j++)
		{
			mpz_mul(prod1, A[row_idx[j]][col_idx[j]], mpz_mat[row_idx[i]][idx_RHS_col]);
			mpz_mul(prod2, A[row_idx[i]][col_idx[j]], mpz_mat[row_idx[j]][idx_RHS_col]);
			mpz_sub(mpz_mat[row_idx[i]][idx_RHS_col], prod1, prod2);

			if (j > 0)
				mpz_divexact(mpz_mat[row_idx[i]][idx_RHS_col], mpz_mat[row_idx[i]][idx_RHS_col], A[row_idx[j - 1]][col_idx[j - 1]]);
		}
	}
}

// Similar as the above function, but the function takes REF-LU matrix, a single vector, and number of rows as arguments; also has option for transpose forward substitution; returns forward-substituted vector
mpz_t* IPGE::FSub(mpz_t **mat, mpz_t *vec, int n, bool transp)
{
	mpz_t* y_temp;
	
	y_temp = GFz_vec_init_zeros(n);
	GFz_vec_copy(vec, y_temp, n);
	
	mpz_t prod1, prod2;
	mpz_inits(prod1, prod2, NULL);
	
	for (int k=0; k < n; k++)
	{	
		for (int i = k+1; i < n; i++)
		{
			if(!transp)
			{
				mpz_mul(prod1, mat[row_idx[k]][col_idx[k]], y_temp[row_idx[i]]);
				mpz_mul(prod2, mat[row_idx[i]][col_idx[k]], y_temp[row_idx[k]]);
				mpz_sub(y_temp[row_idx[i]], prod1, prod2);
				if (k> 0)
					mpz_divexact(y_temp[row_idx[i]], y_temp[row_idx[i]],mat[row_idx[k-1]][col_idx[k-1]]);
			}
			
			else
			{
				mpz_mul(prod1, mat[row_idx[k]][col_idx[k]], y_temp[col_idx[i]]);
				mpz_mul(prod2, mat[row_idx[k]][col_idx[i]], y_temp[col_idx[k]]);
				mpz_sub(y_temp[col_idx[i]], prod1, prod2);
				if (k> 0)
					mpz_divexact(y_temp[col_idx[i]], y_temp[col_idx[i]],mat[row_idx[k-1]][col_idx[k-1]]);
			}
		}	
	}
	return y_temp;
}

/*
	Single FSub step expanding from the 3-argument F-Sub function defined above
	Performs REF Forward Substitution Stepwise Recusion defined in Escobedo (2023) (with the difference that column step_idx of L and previous pivot are not passed as arguments (they are accessed directly)
	The last argument is used to denote whether the FSub is being performed in the transpose sense (i.e., with U^T rather than with L)
*/
void IPGE::FSub_Step(mpz_t **mpz_mat, int mat_cols, int idx_RHS_col, int step_idx, bool transp)
{
	mpz_t prod1, prod2;
	mpz_inits(prod1, prod2, NULL);
	
	
	
	for (int i = step_idx+1; i < rows; i++)
	{
		if(!transp)
		{
			mpz_mul(prod1, A[row_idx[step_idx]][col_idx[step_idx]], mpz_mat[row_idx[i]][idx_RHS_col]);
			mpz_mul(prod2, A[row_idx[i]][col_idx[step_idx]], mpz_mat[row_idx[step_idx]][idx_RHS_col]);
		}
		
		else
		{
			mpz_mul(prod1, A[row_idx[step_idx]][col_idx[step_idx]], mpz_mat[row_idx[i]][idx_RHS_col]);
			mpz_mul(prod2, A[row_idx[step_idx]][col_idx[i]], mpz_mat[row_idx[step_idx]][idx_RHS_col]);
		}
		
		mpz_sub(mpz_mat[row_idx[i]][idx_RHS_col], prod1, prod2);

		if (step_idx > 0)
			mpz_divexact(mpz_mat[row_idx[i]][idx_RHS_col], mpz_mat[row_idx[i]][idx_RHS_col], A[row_idx[step_idx - 1]][col_idx[step_idx - 1]]);
	}
}

/*	
	Single FSub step expanding from the 3-argument F-Sub function defined above
	Performs REF Forward Substitution Stepwise Recusion defined in Escobedo (2023) (with the difference that column step_idx of L and previous pivot are not passed as arguments (they are accessed directly)
	Last argument is used to denote whether the FSub is being performed in the transpose sense (i.e., with U^T rather than with L)
*/
mpz_t* IPGE::FSub_Step_slim(mpz_t *mpz_vec, int step_idx, bool transp)
{
	mpz_t* y_temp;
	
	y_temp = GFz_vec_init_zeros(rows);
	
	mpz_t prod1, prod2;
	mpz_inits(prod1, prod2, NULL);
	
	//Copy the first elements up to step_idx
	for (int i = 0; i <= step_idx; i++)
	{
		if(!transp)
			mpz_set(y_temp[row_idx[i]], mpz_vec[row_idx[i]]);
		else
			mpz_set(y_temp[col_idx[i]], mpz_vec[col_idx[i]]);	
	}
	
	for (int i = step_idx+1; i < rows; i++)
	{
		if(!transp)
		{
			mpz_mul(prod1, A[row_idx[step_idx]][col_idx[step_idx]], mpz_vec[row_idx[i]]);
			mpz_mul(prod2, A[row_idx[i]][col_idx[step_idx]], mpz_vec[row_idx[step_idx]]);
			mpz_sub(y_temp[row_idx[i]], prod1, prod2);
			
			if (step_idx > 0)
				mpz_divexact(y_temp[row_idx[i]], y_temp[row_idx[i]], A[row_idx[step_idx - 1]][col_idx[step_idx - 1]]);
		}
		
		else
		{			
			mpz_mul(prod1, A[row_idx[step_idx]][col_idx[step_idx]],  mpz_vec[col_idx[i]]);
			mpz_mul(prod2, A[row_idx[step_idx]][col_idx[i]],  mpz_vec[col_idx[step_idx]]);
			mpz_sub(y_temp[col_idx[i]], prod1, prod2);
			
			if (step_idx > 0)
			{
				mpz_divexact(y_temp[col_idx[i]], y_temp[col_idx[i]], A[row_idx[step_idx - 1]][col_idx[step_idx - 1]]);
			}
		}		
	}	
	return y_temp;
}

/*
	REF backward substitution, as described in Escobedo and Moreno-Centeno (2015)
	The algorithm solves the scaled triangular system Udet(A)x=det(A)y in exact rational arithmetic 
*/
void IPGE::BSub()
{
	mpz_t prod, det;
	mpz_inits(prod, det, NULL);			

	//Initalize X_p-vector(s) to all 0s
	X_p =  GFz_mat_init_zeros(rows, cols_B);
			
	mpz_set(det, A[row_idx[rows-1]][col_idx[rows-1]]);

	for (int k = 0; k < cols_B; k++)
	{
		for (int i = rows-1; i >= 0; i--)
		{
			//Calculating y':= det(A)y
			mpz_mul(X_p[row_idx[i]][k], det, Y[row_idx[i]][k]);

			for (int j = i+1; j < cols; j++)
			{
				mpz_mul(prod, A[row_idx[i]][col_idx[j]], X_p[row_idx[j]][k]);
				mpz_sub(X_p[row_idx[i]][k], X_p[row_idx[i]][k], prod);
			}

			mpz_divexact(X_p[row_idx[i]][k], X_p[row_idx[i]][k],  A[row_idx[i]][col_idx[i]]);	
		}
	}
}

// Initializes the right-hand side vectors and solution vectors
void IPGE::set_IPGE_Bq_Xq()
{
	//Common denominator to all entries in Xq
	mpq_t det;
	mpq_init(det);
	mpq_set_z(det, A[row_idx[rows-1]][col_idx[rows-1]]);
	
	//Initialize both matrices to 0
	Bq = GFq_mat_init_zeros(rows, cols_B);
	Xq = GFq_mat_init_zeros(rows, cols_B);
	
	//Copy existing B matrix to Bq
	GFzToq_mat_copy(B, Bq, rows, cols_B);

	//Set elements of Xq
	for (int k = 0; k < cols_B; k++)
	{
		for (int i = 0; i < rows; i++)
		{
			mpq_init(Xq[row_idx[i]][k]);
			mpq_set_z(Xq[row_idx[i]][k], X_p[row_idx[i]][k]);
			mpq_div(Xq[row_idx[i]][k], Xq[row_idx[i]][k], det);
		}
	}
}

// Resets all values of the binary frame_sgn array to 1 
void IPGE::reset_frameSigns()
{
	frame_sgn.clear();
	for (int i = 0; i < rows; i++)
		frame_sgn.push_back(1);
}

// REF operation for rewinding all entries within a particular frame (and possibly shifting left)
void IPGE::backtrack(int fr_idx, bool hPart, bool fr_shift)
{
	mpz_t prod1, prod2;
	mpz_inits(prod1, prod2, NULL);

	//Backtrack all entries along the horizontal part of the frame
	if (hPart)
	{
		for (int j = fr_idx; j < col_idx.size(); j++)
		{			
			if (fr_idx > 1)
				mpz_mul(prod1, A[row_idx[fr_idx - 2]][fr_idx - 2], A[row_idx[fr_idx]][col_idx[j]]);
			else
				mpz_set(prod1, A[row_idx[fr_idx]][col_idx[j]]);

			mpz_mul(prod2, A[row_idx[fr_idx - 1]][col_idx[j]], A[row_idx[fr_idx]][col_idx[fr_idx - 1]]);

			mpz_add(prod1, prod1, prod2);

			if (fr_shift)
				mpz_divexact(A[row_idx[fr_idx - 1]][col_idx[j]], prod1, A[row_idx[fr_idx - 1]][col_idx[fr_idx - 1]]);
			else
				mpz_divexact(A[row_idx[fr_idx]][col_idx[j]], prod1, A[row_idx[fr_idx - 1]][col_idx[fr_idx - 1]]);
		}
	}

	//Backtrack all entries along the vertical part of the frame
	else
	{
		for (int i = fr_idx; i < row_idx.size(); i++)
		{
			if (fr_idx > 1)
				mpz_mul(prod1, A[row_idx[fr_idx - 2]][fr_idx - 2], A[row_idx[i]][col_idx[fr_idx]]);
			else
				mpz_set(prod1, A[row_idx[i]][col_idx[fr_idx]]);

			mpz_mul(prod2, A[row_idx[fr_idx - 1]][col_idx[fr_idx]], A[row_idx[i]][col_idx[fr_idx - 1]]);

			mpz_add(prod1, prod1, prod2);

			if (fr_shift)
				mpz_divexact(A[row_idx[i]][col_idx[fr_idx - 1]], prod1, A[row_idx[fr_idx - 1]][col_idx[fr_idx - 1]]);
			else
				mpz_divexact(A[row_idx[i]][col_idx[fr_idx]], prod1, A[row_idx[fr_idx - 1]][col_idx[fr_idx - 1]]);
		}
	}
}

// REF operation for performing a change of pivot operation on a particular frame
void IPGE::RwSOP(int fr_idx, int pivPrime_idx)
{
	mpz_t prod1, prod2;
	mpz_inits(prod1, prod2, NULL);

	for (int j = fr_idx + 1; j < col_idx.size(); j++)
	{
		cout << "[(" << A[row_idx[fr_idx - 1]][col_idx[pivPrime_idx]] << ")(" << A[row_idx[fr_idx]][col_idx[j]] << ")-"
			<< "(" << A[row_idx[fr_idx]][col_idx[pivPrime_idx]] << ")(" << A[row_idx[fr_idx - 1]][col_idx[j]] << ")]"
			<< "/" << A[row_idx[fr_idx - 1]][col_idx[fr_idx - 1]] << "=" << A[row_idx[fr_idx]][col_idx[j]] << endl;

		mpz_mul(prod1, A[row_idx[fr_idx - 1]][col_idx[pivPrime_idx]], A[row_idx[fr_idx]][col_idx[j]]);
		mpz_mul(prod2, A[row_idx[fr_idx]][col_idx[pivPrime_idx]], A[row_idx[fr_idx - 1]][col_idx[j]]);
		mpz_sub(prod1, prod1, prod2);
		mpz_divexact(A[row_idx[fr_idx]][col_idx[j]], prod1, A[row_idx[fr_idx - 1]][col_idx[fr_idx - 1]]);
	}
}

// Backtrack subroutine specific to column-push operation
void IPGE::CPush_backtrack(int fr_idx)
{
	mpz_t prod1, prod2;
	mpz_inits(prod1, prod2, NULL);

	//Backtrack all entries along the vertical part of the frame and shift to the adjacent frame to the left
	
	for (int i = fr_idx; i < row_idx.size(); i++)
	{
		if (fr_idx > 1)
			mpz_mul(prod1, A[row_idx[fr_idx - 2]][col_idx[fr_idx - 2]], A[row_idx[i]][col_idx[fr_idx]]);
		else
			mpz_set(prod1, A[row_idx[i]][col_idx[fr_idx]]);

		mpz_mul(prod2, A[row_idx[fr_idx - 1]][col_idx[fr_idx]], A[row_idx[i]][col_idx[fr_idx - 1]]);

		mpz_add(prod1, prod1, prod2);

		mpz_divexact(A[row_idx[i]][col_idx[fr_idx-1]], prod1, A[row_idx[fr_idx - 1]][col_idx[fr_idx - 1]]);
	}
}

// Backtrack subroutine specific to column-push operation
void IPGE::CPush_backtrack_slim(int fr_idx)
{
	mpz_t prod1, prod2;
	mpz_inits(prod1, prod2, NULL);

	//Backtrack all entries along the vertical part of the frame and shift to the adjacent frame to the left

	for (int i = fr_idx; i < row_idx.size(); i++)
	{
		if (fr_idx > 1)
			mpz_mul(prod1, A[row_idx[fr_idx - 2]][col_idx[fr_idx - 2]], A[row_idx[i]][col_idx[fr_idx]]);
		else
			mpz_set(prod1, A[row_idx[i]][col_idx[fr_idx]]);

		//When current entry being backtracked should have the opposite sign of the one currently stored
		if (frame_sgn[fr_idx] == -1)
			mpz_neg(prod1, prod1);

		mpz_mul(prod2, A[row_idx[fr_idx - 1]][col_idx[fr_idx]], A[row_idx[i]][col_idx[fr_idx - 1]]);

		mpz_add(prod1, prod1, prod2);

		mpz_divexact(A[row_idx[i]][col_idx[fr_idx - 1]], prod1, A[row_idx[fr_idx - 1]][col_idx[fr_idx - 1]]);
	}
}

// RwSOP subroutine specific to column-push operation
void IPGE::CPush_RwSOP(int fr_idx, int c_plus_idx)
{
	mpz_t prod1, prod2;
	mpz_inits(prod1, prod2, NULL);

	for (int j = fr_idx+1; j < col_idx.size(); j++)
	{		
		mpz_mul(prod1, A[row_idx[fr_idx - 1]][col_idx[fr_idx]], A[row_idx[fr_idx]][col_idx[j]]);
		mpz_mul(prod2, A[row_idx[fr_idx]][col_idx[fr_idx]], A[row_idx[fr_idx-1]][col_idx[j]]);
		mpz_sub(prod1, prod1, prod2);
		mpz_divexact(A[row_idx[fr_idx]][col_idx[j]], prod1, A[row_idx[fr_idx - 1]][col_idx[fr_idx - 1]]);
	}

	//Apply operations to entering column
	mpz_mul(prod1, A[row_idx[fr_idx - 1]][col_idx[fr_idx]], A_plus[row_idx[fr_idx]][c_plus_idx]);
	mpz_mul(prod2, A[row_idx[fr_idx]][col_idx[fr_idx]], A_plus[row_idx[fr_idx - 1]][c_plus_idx]);
	mpz_sub(prod1, prod1, prod2);
	mpz_divexact(A_plus[row_idx[fr_idx]][c_plus_idx], prod1, A[row_idx[fr_idx - 1]][col_idx[fr_idx - 1]]);
}

// Reduced RwSOP subroutine specific to column-push operation
void IPGE::CPush_RwSOP_slim(int fr_idx, int c_plus_idx)
{
	mpz_t prod1, prod2;
	mpz_inits(prod1, prod2, NULL);

	for (int j = fr_idx + 1; j < col_idx.size(); j++)
	{
		mpz_mul(prod1, A[row_idx[fr_idx - 1]][col_idx[fr_idx]], A[row_idx[fr_idx]][col_idx[j]]);
		mpz_mul(prod2, A[row_idx[fr_idx]][col_idx[fr_idx]], A[row_idx[fr_idx - 1]][col_idx[j]]);
		mpz_sub(prod1, prod1, prod2);
		mpz_divexact(A[row_idx[fr_idx]][col_idx[j]], prod1, A[row_idx[fr_idx - 1]][col_idx[fr_idx - 1]]);

		//When current entry being backtracked should have the opposite sign of the one currently stored
		if (frame_sgn[fr_idx] == -1)
			mpz_neg(A[row_idx[fr_idx]][col_idx[j]], A[row_idx[fr_idx]][col_idx[j]]);
	}

	//Apply operations to entering column
	mpz_mul(prod1, A[row_idx[fr_idx - 1]][col_idx[fr_idx]], A_plus[row_idx[fr_idx]][c_plus_idx]);
	mpz_mul(prod2, A[row_idx[fr_idx]][col_idx[fr_idx]], A_plus[row_idx[fr_idx - 1]][c_plus_idx]);
	mpz_sub(prod1, prod1, prod2);
	mpz_divexact(A_plus[row_idx[fr_idx]][c_plus_idx], prod1, A[row_idx[fr_idx - 1]][col_idx[fr_idx - 1]]);

	//When current entry being backtracked should have the opposite sign of the one currently stored
	if (frame_sgn[fr_idx] == -1)
		mpz_neg(A_plus[row_idx[fr_idx]][c_plus_idx], A_plus[row_idx[fr_idx]][c_plus_idx]);
}

// Column permutation subroutine specific to column-push operation
void IPGE::CPush_CPerm(int fr_idx)
{
	mpz_t temp;
	mpz_init(temp);

	//For frames 1,...,fr_idx-1, perform a normal permutation of columns fr_idx and fr_idx-1
	std::swap(col_idx[fr_idx], col_idx[fr_idx-1]);

	//For the remaining frames, change sign of entries along frame fr_idx-1 (fr_idx before abvoe column permuation) and then re-permute the entries along columns fr_idx and fr_idx-1
	for (int i = fr_idx; i < rows; i++)
	{
		mpz_neg(temp, A[row_idx[i]][col_idx[fr_idx - 1]]);
		mpz_set(A[row_idx[i]][col_idx[fr_idx - 1]], A[row_idx[i]][col_idx[fr_idx]]);
		mpz_set(A[row_idx[i]][col_idx[fr_idx]], temp);
	}
}

// Reduced column permutation subroutine specific to column-push operation
void IPGE::CPush_CPerm_slim(int fr_idx)
{
	mpz_t temp;
	mpz_init(temp);

	//For frames 1,...,fr_idx-1, perform a normal permutation of columns fr_idx and fr_idx-1
	std::swap(col_idx[fr_idx], col_idx[fr_idx - 1]);

	//For the remaining frames, change sign of entries along frame fr_idx-1 (fr_idx before abvoe column permuation) and then re-permute the entries along columns fr_idx and fr_idx-1
	for (int i = fr_idx; i < rows; i++)
	{
		//When current frame matches the sign of the stored entry
		if (frame_sgn[fr_idx] == 1)
			mpz_neg(temp, A[row_idx[i]][col_idx[fr_idx - 1]]);

		//When current frame does not match the sign of the stored entry
		else
			mpz_set(temp, A[row_idx[i]][col_idx[fr_idx - 1]]);

		mpz_set(A[row_idx[i]][col_idx[fr_idx - 1]], A[row_idx[i]][col_idx[fr_idx]]);
		mpz_set(A[row_idx[i]][col_idx[fr_idx]], temp);
	}
}

// Changes signs of ALL entries in frames k,...,n (includes extension to entering column)
void IPGE::CPush_flipSigns(int fr_idx, int c_plus_idx)
{
	for (int i = fr_idx; i < rows; i++)
	{
		for (int j = fr_idx; j < cols; j++)
			mpz_neg(A[row_idx[i]][col_idx[j]], A[row_idx[i]][col_idx[j]]);
		
		//Apply operations to entering column
		mpz_neg(A_plus[row_idx[i]][c_plus_idx], A_plus[row_idx[i]][c_plus_idx]);
	}	
}

// Reduced function for changing signs of ALL entries in frames k,...,n (includes extension to entering column)
void IPGE::CPush_flipSigns_slim(int fr_idx)
{
	GF_vec_consencutive_sgn_flip(frame_sgn, rows, fr_idx);
}

// Column-push algorithm 
void IPGE::CPush(int col_idx_right)
{
	int c_plus_idx = ex_idx - col_idx_right + num_ent;

	for (int i = col_idx_right; i < rows; i++)
	{
		CPush_backtrack(i);
		CPush_RwSOP(i, c_plus_idx);
		CPush_CPerm(i);
		if (i+1 < rows)
			CPush_flipSigns(i+1, c_plus_idx);
	}

	//Swap exiting column with entering column
	for (int i = 0; i < rows; i++)
		mpz_set(A[row_idx[i]][col_idx[rows-1]], A_plus[row_idx[i]][c_plus_idx]);
}

// Reduced column-push algorithm 
void IPGE::CPush_slim(int col_idx_right)
{
	reset_frameSigns();
	
	int c_plus_idx = ex_idx - col_idx_right + num_ent;

	//cout << "c_plus_idx:" << c_plus_idx + 1 << endl;

	for (int i = col_idx_right; i < rows; i++)
	{
		CPush_backtrack_slim(i);
		CPush_RwSOP_slim(i, c_plus_idx);
		CPush_CPerm_slim(i);
		CPush_flipSigns_slim(i + 1);
	}

	//Swap exiting column with entering column
	for (int i = 0; i < rows; i++)
		mpz_set(A[row_idx[i]][col_idx[rows - 1]], A_plus[row_idx[i]][c_plus_idx]);
}

// Reduced column-push algorithm 
void IPGE::CPush_slim(int col_idx_right, bool printIt)
{
	reset_frameSigns();
	
	int c_plus_idx = ex_idx - col_idx_right + num_ent;

	for (int i = col_idx_right; i < rows; i++)
	{
		CPush_backtrack_slim(i);
		CPush_RwSOP_slim(i, c_plus_idx);
		CPush_CPerm_slim(i);
		CPush_flipSigns_slim(i + 1);
		
		if(printIt)
			GFz_idx_formPrint(A, row_idx, col_idx, pSpace + 3, 0);
	}

	//Swap exiting column with entering column
	for (int i = 0; i < rows; i++)
		mpz_set(A[row_idx[i]][col_idx[rows - 1]], A_plus[row_idx[i]][c_plus_idx]);
}

// Push-and-swap REF column replacement update
void IPGE::Update()
{
	A_plus = GFz_mat_init_density(rows, num_ent, get_lb(), get_ub(), get_seed1() - get_seed2() * 5, 4 * get_seed1() + get_seed2(), get_density(), false);

	vector<int> temp;
	for (int i = 0; i < num_ent; i++)
		temp.push_back(i);

	for (int i = ex_idx + num_ent; i >= ex_idx + 1; i--)
	{
		FSub(A_plus, num_ent, num_ent - (i - ex_idx));
		CPush(i);		
	}
}

// Reduced push-and-swap REF column replacement update
void IPGE::Update_slim()
{
	num_ent = 1;

	A_plus = GFz_mat_init_density(rows, num_ent, get_lb(), get_ub(), get_seed1() - get_seed2() * 5, 4 * get_seed1() + get_seed2(), get_density(), false);
	A_plus0 = GFz_mat_init_density(rows, num_ent, get_lb(), get_ub(), get_seed1() - get_seed2() * 5, 4 * get_seed1() + get_seed2(), get_density(), false);


	vector<int> temp;
	for (int i = 0; i < num_ent; i++)
		temp.push_back(i);

	for (int i = ex_idx + num_ent; i >= ex_idx + 1; i--)
	{
		reset_frameSigns();
		FSub(A_plus, num_ent, num_ent - (i - ex_idx));
		CPush_slim(i);
	}
}



// Changes signs of ALL entries in frames k,...,n, for rank-one update
void IPGE::ROU_flipSigns(int fr_idx)
{
	for (int i = fr_idx; i < rows; i++)
	{
		for (int j = fr_idx; j < cols; j++)
			mpz_neg(A[row_idx[i]][col_idx[j]], A[row_idx[i]][col_idx[j]]);
	}	
}

// Changes signs of rank-one update vector starting at the specified index
void IPGE::ROU_flipSigns_vec(mpz_t * mpz_vec, int beg_idx)
{
	for (int i = beg_idx; i < rows; i++)
		mpz_neg(mpz_vec[row_idx[i]], mpz_vec[row_idx[i]]);
}

// RwSOP subroutine specific to rank-one update
void IPGE::ROU_RwSOP(int fr_idx)
{
	mpz_t prod1, prod2;
	mpz_inits(prod1, prod2, NULL);

	for (int j = fr_idx + 1; j < col_idx.size(); j++)
	{
		mpz_mul(prod1, A[row_idx[fr_idx - 1]][col_idx[fr_idx]], A[row_idx[fr_idx]][col_idx[j]]);
		mpz_mul(prod2, A[row_idx[fr_idx]][col_idx[fr_idx]], A[row_idx[fr_idx - 1]][col_idx[j]]);
		mpz_sub(prod1, prod1, prod2);
		mpz_divexact(A[row_idx[fr_idx]][col_idx[j]], prod1, A[row_idx[fr_idx - 1]][col_idx[fr_idx - 1]]);
	}
}

// Backtrack subroutine specific to rank-one update
void IPGE::ROU_backtrack(int fr_idx)
{
	mpz_t prod1, prod2;
	mpz_inits(prod1, prod2, NULL);

	//Backtrack all entries along the vertical part of the frame and shift to the adjacent frame to the left
	for (int i = fr_idx; i < row_idx.size(); i++)
	{
		if (fr_idx > 1)
			mpz_mul(prod1, A[row_idx[fr_idx - 2]][col_idx[fr_idx - 2]], A[row_idx[i]][col_idx[fr_idx]]);
		else
			mpz_set(prod1, A[row_idx[i]][col_idx[fr_idx]]);


		mpz_mul(prod2, A[row_idx[fr_idx - 1]][col_idx[fr_idx]], A[row_idx[i]][col_idx[fr_idx - 1]]);
		mpz_add(prod1, prod1, prod2);
		mpz_divexact(A[row_idx[i]][col_idx[fr_idx - 1]], prod1, A[row_idx[fr_idx - 1]][col_idx[fr_idx - 1]]);
	}
}

//Sets the update vector v either to either have all entries randomly generated (Experiment 1 in the paper) or as partially as a linear combination of sub-columns of A0[i][j] (Experiment 2 in the paper)
mpz_t* IPGE::ROU_set_update_vec_v(int dep_col_idx, int dep_row_idx_end, bool isLinearCombination)
{
	mpz_t temp, ZERO;
	
	//Initialize a vector with random entries, which will be used to generate some of the entries of the update vector v (those that are not linearly dependent)
	mpz_t *v_temp, *v_temp_permuted;
	v_temp = GFz_vec_init_zeros(rows);
	v_temp = GFz_vec_init_rand(rows, get_lb(), get_ub(), get_seed1() - get_seed2() * 5);
	
	mpz_init(ZERO);
	
	//If indicated by the isLinearCombination argument, initialize the rank-one update vector v as a linear combination of the columns of A0
	if(isLinearCombination)
	{		
		//Reset first dep_idx_end entries to zero since these will be equal to one of the input columns
		for(int i=0; i <= dep_row_idx_end; i++)
			mpz_set(v_temp[i], ZERO);
		
		//This for-loop either adds (when j==dep_col_idx) or substracts (when j!=dep_col_idx) the subvector consisting of rows 0...dep_idx_end of column j from the original input matrix (A0[i][j]) 
		for (int j=0; j < dep_col_idx; j++)
		{
			for(int i=0; i <= dep_row_idx_end; i++)
			{				
				if(j == dep_col_idx)
					mpz_add(v_temp[i], v_temp[i], A0[i][j]);
				else				
					mpz_sub(v_temp[i], v_temp[i], A0[i][j]);				
			}						
		}		
	}
	
	//Initialize v_temp_permuted according to the row/column permutations of the current REF LU factorization
	v_temp_permuted = GFz_vec_init_zeros(rows);
	
	for (int i=0; i < rows; i++)
		mpz_set(v_temp_permuted[row_idx[i]], v_temp[i]);
 
	return v_temp_permuted;
}

//Sets the update vector w to have all entries randomly generated (for Experiments 1 and 2 in the paper) 
mpz_t* IPGE::ROU_set_update_vec_w()
{
	mpz_t *w_temp, *w_temp_permuted;
	
	//Initialize rank-one update vector randomly	
	w_temp = GFz_vec_init_rand(cols, get_lb(), get_ub(), get_seed1() + 2*get_seed2());
	
	w_temp_permuted = GFz_vec_init_zeros(cols);
	
	for (int j=0; j < cols; j++)
		mpz_set(w_temp_permuted[col_idx[j]],w_temp[j]);
	
	return w_temp_permuted;
}

// Checks that the next divisor needed to obtain entries of the updated REF-LU factorization is nonzero
bool IPGE::ROU_checkNextDivisorGood(mpz_t **mat, mpz_t* mpz_vec, int k, bool transp)
{
	mpz_t prod1, prod2, temp;
	mpz_inits(prod1, prod2, NULL);
	mpz_init(temp);
	
	bool good_divisor = true;
		
	if(!transp)
	{
		mpz_mul(prod1, mat[row_idx[k-1]][col_idx[k-1]], mpz_vec[row_idx[k]]);
		mpz_mul(prod2, mat[row_idx[k]][col_idx[k-1]], mpz_vec[row_idx[k-1]]);	
	}
	
	else
	{
		mpz_mul(prod1, mat[row_idx[k-1]][col_idx[k-1]], mpz_vec[col_idx[k]]);
		mpz_mul(prod2, mat[row_idx[k-1]][col_idx[k]], mpz_vec[col_idx[k-1]]);
	}
	
	if (mpz_cmp(prod1, prod2)==0)
			good_divisor = false;
	
	return good_divisor;
}

// Recalculates an entry from REF-LU(Ah0) up to the specified iteration k
void IPGE::ROU_recomputeAhDiagEntry(mpz_t **mpz_mat, int i, int j, int k)
{
	mpz_t prod1, prod2;
	mpz_inits(prod1, prod2, NULL);
	
	for (int step=0; step < k-1; step++)
	{
		mpz_mul(prod1, Ah[row_idx[step]][col_idx[step]], 
								Ah[row_idx[i]][col_idx[j]]);

		mpz_mul(prod2, Ah[row_idx[step]][col_idx[j]], 
						Ah[row_idx[i]][col_idx[step]]);

		mpz_sub(Ah[row_idx[i]][col_idx[j]], prod1, prod2);
		
		if(step > 0)			
			mpz_divexact (Ah[row_idx[i]][col_idx[j]], Ah[row_idx[i]][col_idx[j]], Ah[row_idx[step-1]][col_idx[step-1]]);
	}	
}

// Handles the case when the divisor needed to obtain entries of L is zero; it is a subroutine called, as needed, within the ROU algorithm
void IPGE::ROU_Special_Case(int k, bool printIt)
{
	//Performs REF adjacent column permutation on REF-LU(A)
	mpz_t temp, ZERO, prod1,prod2;
			
	mpz_inits(temp,ZERO, NULL);
	mpz_inits(prod1,prod2,NULL);
	
	
	//Perform adjacent column permutation
	ROU_backtrack(k);
	
	if(printIt)
	{
		cout << "After backtracking:" << endl;
		GFz_idx_formPrint(A, row_idx, col_idx, pSpace + 3, 0);
	}
	
	ROU_RwSOP(k);
	if(printIt)
	{
		cout << "After RwSop:" << endl;
		GFz_idx_formPrint(A, row_idx, col_idx, pSpace + 3, 0);
	}

	CPush_CPerm(k);
	if(printIt)
	{
		cout << "After CPerm:" << endl;
		GFz_idx_formPrint(A, row_idx, col_idx, pSpace + 3, 0);
	}

	ROU_flipSigns(k+1);
	
	if(printIt)
	{
		cout << "After Sign flips:" << endl; 
		GFz_idx_formPrint(A, row_idx, col_idx, pSpace + 3, 0);	
	}
	
	//Places correct entry in place k,k
	if(k==2)
	{
		mpz_init(temp);
		mpz_set(temp, A0[row_idx[k]][col_idx[k]]);
		mpz_mul(prod1, v[row_idx[k]], w[col_idx[k]]);
		mpz_add(Ah[row_idx[k]][col_idx[k]],temp,prod1); 
	}
	
	if (k >2)
	{
		mpz_init(temp);
		mpz_set(temp, A0[row_idx[k]][col_idx[k]]);
		mpz_mul(prod1, v[row_idx[k]], w[col_idx[k]]);
		mpz_add(temp, temp, prod1);	
		
		mpz_set(Ah[row_idx[k]][col_idx[k]],temp); 	
		
		ROU_recomputeAhDiagEntry(Ah, k, k, k-1);
	}
		
	if(k >=2)
	{
		//Calculate k-1,k-1 entry
		mpz_init(temp);
		mpz_set(temp, A0[row_idx[k-1]][col_idx[k-1]]);
		mpz_mul(prod1, v[row_idx[k-1]], w[col_idx[k-1]]);
		mpz_add(temp, temp, prod1);	
		
		if(printIt)
		{
			cout << "Initial entry = " << A0[row_idx[k-1]][col_idx[k-1]] << " + ";
			cout << "(" << v[row_idx[k-1]] << ")(" << w[col_idx[k-1]] << ") = " << temp << "\n";
		}

		mpz_set(Ah[row_idx[k-1]][col_idx[k-1]],temp); 		
		
		if(printIt)
			GFz_idx_formPrint(Ah, row_idx, col_idx, pSpace, 0);
		
		ROU_recomputeAhDiagEntry(Ah, k-1, k-1, k-1);
	}
}

/*
	This is the main algorithm discussed in the paper, with adjustments for special case 2 (SC2) included
	This function returns an int array consisting of (1) the column index used for generating v, (2) the row index used for generating v, and (3) the number of SC2 calls required
*/
int* IPGE::Rank_One_Update(bool printIt, int up_idx_1, int up_idx_2)
{
	//Needed for calculations
	mpz_t temp, ZERO, prod1,prod2;
			
	mpz_inits(temp,ZERO, NULL);
	mpz_inits(prod1,prod2,NULL);
	
	//Values to be returned by this function
	static int rou_outputs [3];
	
	//Keep track of times columns are permuted during the algorithm to avoid dividing by zero
	num_perms = 0;
	
	//The update vector w is initialized to have all randomized entries (within the specified integer ranges), that is, it is not generated from a linear combination of the input matrix
	w = ROU_set_update_vec_w();
	
	//This indicates that update vector v should also be generated to have all randomized entries (within the specified integer ranges)
	if(up_idx_1 < -1)
	{
		v = ROU_set_update_vec_v(up_idx_1, up_idx_2, 0); //Since v is not designed to be linearly dependent on the columns of A, give 0 as the last argument
		rou_outputs[0] = -2;
		rou_outputs[1] = -2;
	}

	//This indicates that v should be set to force its first up_idx_2 elements to be a linear combination of the first up_idx_1 columns of the input matrix
	else
	{
		if(up_idx_1 == -1)//Sets up_idx_1 and up_idx_2 randomly (as is done in Experiment 2 of Escobedo (2013)
		{
			srand(-4*get_seed1()+cols*get_seed2());

			int idx_random = rand() % cols;    
						
			up_idx_1 = idx_random; 
			up_idx_2 = rand() % (cols-idx_random) + idx_random; 	
		} 
		v = ROU_set_update_vec_v(up_idx_1, up_idx_2, 1); //Since v is designed to be linearly dependent on the columns of A, give 1 as the last argument
		
		rou_outputs[0] = up_idx_1; 
		rou_outputs[1] = up_idx_2;
	}
	
	if(printIt)
	{
		cout << "Update vector v:\n";
		GF_idx_rawPrint(v, row_idx);
		cout << "Update vector w:\n";
		GF_idx_rawPrint(w, col_idx);
	}

	//Initialize the four forward sub vectors needed by the algorithm
	y_k = GFz_vec_init_zeros(rows);
	y_k_old = GFz_vec_init_zeros(rows);
	z_k = GFz_vec_init_zeros(cols);
	z_k_old = GFz_vec_init_zeros(cols);

	// Set y and z to v and w, respectively
	for (int i = 0; i < rows; i++)
		mpz_set(y_k_old[row_idx[i]], v[row_idx[i]]);

	for (int j = 0; j < cols; j++)
		mpz_set(z_k_old[col_idx[j]], w[col_idx[j]]);
	
	// Check if element that will be used as divisor to get second column of L is zero
	if(ROU_checkNextDivisorGood(A, y_k_old, 1, false)==0)
	{
		if(printIt)
			cout << "ZERO DIVISOR detected\n";
		ROU_Special_Case(1,printIt);
		num_perms++;
	}
	// Perform first iteration of REF forward sub
	y_k = FSub_Step_slim(y_k_old, 0, false);
	z_k = FSub_Step_slim(z_k_old, 0, true);

	//Initalize matrix that will hold the updated REF-LU factorization 
	Ah = GFz_mat_init_zeros(rows, cols);

	//Obtain first column of Ah; also, obtain original diagonal of Ah
	for (int i = 0; i < rows; i++)
	{
		mpz_mul(temp, v[row_idx[i]], w[col_idx[i]]);
		mpz_add(Ah[row_idx[i]][col_idx[i]], A0[row_idx[i]][col_idx[i]],temp);
		
		if(i>0)
		{
			mpz_mul(temp, v[row_idx[i]], w[col_idx[0]]);
			mpz_add(Ah[row_idx[i]][col_idx[0]],A[row_idx[i]][col_idx[0]],temp); 
		}
	}

	//Obtain first row of Ah; 
	for (int j = 1; j < cols; j++)
	{
		mpz_mul(temp, v[row_idx[0]], w[col_idx[j]]);
		mpz_add(Ah[row_idx[0]][col_idx[j]],A[row_idx[0]][col_idx[j]],temp); 	
	}

	
	if(printIt)
	{
		cout << "ROU: Output from initial steps\n";
		GFz_idx_formPrint(Ah, row_idx, col_idx, pSpace, 0);
	}
	
	for (int k = 1; k < rows-1; k++)
	{
		if(printIt)
			cout << "Inner for-loop iteration " << k+1 << "\n";
		GFz_vec_copy(y_k, y_k_old, rows);	
		GFz_vec_copy(z_k, z_k_old, cols);
		
		if(ROU_checkNextDivisorGood(A, y_k_old, k+1, false)==0)
		{
			if(printIt)
				cout << "ZERO DIVISOR detected\n";
			ROU_Special_Case(k+1,printIt);
			num_perms++;
		}

		if(printIt)
			GFz_idx_formPrint(Ah, row_idx, col_idx, pSpace, 0);
		
		//Perform next REF FS iteration
		y_k = FSub_Step_slim(y_k, k, false);
		z_k = FSub_Step_slim(z_k, k, true);
						
		for (int i = k; i < rows; i++)
		{
			//Calculate next IPGE iteration of Ah diagonals
			mpz_mul(prod1, Ah[row_idx[k-1]][col_idx[k-1]], Ah[row_idx[i]][col_idx[i]]);
			mpz_mul(prod2, Ah[row_idx[k-1]][col_idx[i]], Ah[row_idx[i]][col_idx[k-1]]);
			
			mpz_sub(Ah[row_idx[i]][col_idx[i]], prod1,prod2);
			
			if (k>1)
				mpz_divexact(Ah[row_idx[i]][col_idx[i]],Ah[row_idx[i]][col_idx[i]], Ah[row_idx[k-2]][col_idx[k-2]]);
		}
		
		//Get next column of Lh matrix
		for (int i = k+1; i < rows; i++)
		{			
			mpz_mul(prod1, Ah[row_idx[k]][col_idx[k]], y_k_old[row_idx[i]]);
			mpz_mul(prod2, Ah[row_idx[k-1]][col_idx[k-1]], y_k[row_idx[i]]);
			mpz_sub(Ah[row_idx[i]][col_idx[k]], prod1,prod2);
			
			mpz_divexact(Ah[row_idx[i]][col_idx[k]],Ah[row_idx[i]][col_idx[k]], y_k_old[row_idx[k]]);
		}
		
		//Get next row of Uh matrix
		for (int j = k+1; j < cols; j++)
		{
			mpz_mul(prod1, Ah[row_idx[k]][col_idx[k]], z_k_old[col_idx[j]]);
			mpz_mul(prod2, Ah[row_idx[k-1]][col_idx[k-1]], z_k[col_idx[j]]);
			mpz_sub(Ah[row_idx[k]][col_idx[j]], prod1,prod2);
			
			mpz_divexact(Ah[row_idx[k]][col_idx[j]],Ah[row_idx[k]][col_idx[j]], z_k_old[col_idx[k]]);
		}
		
		if(printIt)
			GFz_idx_formPrint(Ah, row_idx, col_idx, pSpace, 0);
	} 	
	
	//Get last updated diagonal element of Ah
	mpz_mul(prod1, Ah[row_idx[rows-2]][col_idx[rows-2]], Ah[row_idx[rows-1]][col_idx[rows-1]]);
	mpz_mul(prod2, Ah[row_idx[rows-2]][col_idx[rows-1]],Ah[row_idx[rows-1]][col_idx[rows-2]]);
	mpz_sub(Ah[row_idx[rows-1]][col_idx[rows-1]], prod1,prod2);
	
	mpz_divexact(Ah[row_idx[rows-1]][col_idx[rows-1]],Ah[row_idx[rows-1]][col_idx[rows-1]], Ah[row_idx[rows-3]][col_idx[rows-3]]);
	
	if(printIt)
	{
		cout << "Final factorization:\n";
		GFz_idx_formPrint(Ah, row_idx, col_idx, pSpace, 0);
	}
	
	rou_outputs[2] = num_perms;
	
	return rou_outputs;
}

/*
	This is the special variant of the ROU algorithm for performing a column replacement update discussed in the paper, with adjustments for special case 2 (SC2) included 
	The number of SC2 calls required by the algorithm is returned by this function
*/
int* IPGE::Rank_One_Update_CR(bool printIt, int ex_idx)
{
	//Needed for calculations
	mpz_t temp, ZERO, prod1,prod2;
			
	mpz_inits(temp,ZERO, NULL);
	mpz_inits(prod1,prod2,NULL);
	
	
	//Values to be returned by this function; for this special use of the ROU algorithm, the first 
	static int rou_outputs [3];
	
	//The first two elements are set to -2 since special initialization of v is not relevant
	rou_outputs[0] = -2;
	rou_outputs[1] = -2;
	
	
	// Third element keeps track of times columns are permuted to avoid dividing by zero (i.e., SC2 calls)
	rou_outputs[2] = 0;
			
	//Initialize rank-one update vectors	
	
	v = GFz_vec_init_zeros(rows);
	for (int i = 0; i < rows; i++)
		mpz_set(v[row_idx[i]], A_plus0[row_idx[i]][ex_idx]);
	
	w = GFz_vec_init_zeros(cols);
	mpz_set_ui(w[ex_idx],1); 
	
	if(printIt)
	{
		cout << "Update vector v:\n";
		GF_idx_rawPrint(v, row_idx);
		cout << "Update vector w:\n";
		GF_idx_rawPrint(w, col_idx);
	}

	//Initialize the four forward sub vectors needed by the algorithm
	y_k = GFz_vec_init_zeros(rows);
	y_k_old = GFz_vec_init_zeros(rows);
	z_k = GFz_vec_init_zeros(cols);
	z_k_old = GFz_vec_init_zeros(cols);

	// Set y and z to v and w, respectively
	for (int i = 0; i < rows; i++)
		mpz_set(y_k_old[row_idx[i]], v[row_idx[i]]);

	for (int j = 0; j < cols; j++)
		mpz_set(z_k_old[col_idx[j]], w[col_idx[j]]);
	
	// Check if element that will be used as divisor to get second column of L is zero
	if(ROU_checkNextDivisorGood(A, y_k_old, 1, false)==0)
	{
		if(printIt)
			cout << "ZERO DIVISOR detected\n";
		ROU_Special_Case(1,printIt);
		rou_outputs[2]++;
	}
	// Perform first iteration of REF forward sub
	y_k = FSub_Step_slim(y_k_old, 0, false);
	z_k = FSub_Step_slim(z_k_old, 0, true);

	//Initalize matrix that will hold the updated REF-LU factorization 
	Ah = GFz_mat_init_zeros(rows, cols);

	//Obtain first column of Ah; also, obtain original diagonal of Ah
	for (int i = 0; i < rows; i++)
	{
		mpz_mul(temp, v[row_idx[i]], w[col_idx[i]]);
		mpz_add(Ah[row_idx[i]][col_idx[i]], A0[row_idx[i]][col_idx[i]],temp);
		
		if(i>0)
		{
			mpz_mul(temp, v[row_idx[i]], w[col_idx[0]]);
			mpz_add(Ah[row_idx[i]][col_idx[0]],A[row_idx[i]][col_idx[0]],temp); 
		}
	}

	//Obtain first row of Ah; 
	for (int j = 1; j < cols; j++)
	{
		mpz_mul(temp, v[row_idx[0]], w[col_idx[j]]);
		mpz_add(Ah[row_idx[0]][col_idx[j]],A[row_idx[0]][col_idx[j]],temp); 	
	}

	
	if(printIt)
	{
		cout << "ROU: Output from initial steps\n";
		GFz_idx_formPrint(Ah, row_idx, col_idx, pSpace, 0);
	}
	
	for (int k = 1; k < rows-1; k++)
	{
		if(printIt)
			cout << "Inner for-loop iteration " << k+1 << "\n";
		GFz_vec_copy(y_k, y_k_old, rows);	
		GFz_vec_copy(z_k, z_k_old, cols);
		
		if(ROU_checkNextDivisorGood(A, y_k_old, k+1, false)==0)
		{
			if(printIt)
				cout << "ZERO DIVISOR detected\n";
			ROU_Special_Case(k+1,printIt);
			rou_outputs[2]++;
		}

		if(printIt)
			GFz_idx_formPrint(Ah, row_idx, col_idx, pSpace, 0);
		
		//Perform next REF FS iteration
		y_k = FSub_Step_slim(y_k, k, false);
		z_k = FSub_Step_slim(z_k, k, true);
						
		for (int i = k; i < rows; i++)
		{
			//Calculate next IPGE iteration of Ah diagonals
			mpz_mul(prod1, Ah[row_idx[k-1]][col_idx[k-1]], Ah[row_idx[i]][col_idx[i]]);
			mpz_mul(prod2, Ah[row_idx[k-1]][col_idx[i]], Ah[row_idx[i]][col_idx[k-1]]);
			
			mpz_sub(Ah[row_idx[i]][col_idx[i]], prod1,prod2);
			
			if (k>1)
				mpz_divexact(Ah[row_idx[i]][col_idx[i]],Ah[row_idx[i]][col_idx[i]], Ah[row_idx[k-2]][col_idx[k-2]]);
		}
		
		//Get next column of Lh matrix
		for (int i = k+1; i < rows; i++)
		{			
			mpz_mul(prod1, Ah[row_idx[k]][col_idx[k]], y_k_old[row_idx[i]]);
			mpz_mul(prod2, Ah[row_idx[k-1]][col_idx[k-1]], y_k[row_idx[i]]);
			mpz_sub(Ah[row_idx[i]][col_idx[k]], prod1,prod2);
			
			mpz_divexact(Ah[row_idx[i]][col_idx[k]],Ah[row_idx[i]][col_idx[k]], y_k_old[row_idx[k]]);
		}
		
		//Get next row of Uh matrix
		for (int j = k+1; j < cols; j++)
		{
			mpz_mul(prod1, Ah[row_idx[k]][col_idx[k]], z_k_old[col_idx[j]]);
			mpz_mul(prod2, Ah[row_idx[k-1]][col_idx[k-1]], z_k[col_idx[j]]);
			mpz_sub(Ah[row_idx[k]][col_idx[j]], prod1,prod2);
			
			mpz_divexact(Ah[row_idx[k]][col_idx[j]],Ah[row_idx[k]][col_idx[j]], z_k_old[col_idx[k]]);
		}
		
		if(printIt)
			GFz_idx_formPrint(Ah, row_idx, col_idx, pSpace, 0);
	} 	
	
	//Get last updated diagonal element of Ah
	mpz_mul(prod1, Ah[row_idx[rows-2]][col_idx[rows-2]], Ah[row_idx[rows-1]][col_idx[rows-1]]);
	mpz_mul(prod2, Ah[row_idx[rows-2]][col_idx[rows-1]],Ah[row_idx[rows-1]][col_idx[rows-2]]);
	mpz_sub(Ah[row_idx[rows-1]][col_idx[rows-1]], prod1,prod2);
	
	mpz_divexact(Ah[row_idx[rows-1]][col_idx[rows-1]],Ah[row_idx[rows-1]][col_idx[rows-1]], Ah[row_idx[rows-3]][col_idx[rows-3]]);
	
	if(printIt)
	{
		cout << "Final factorization:\n";	
		GFz_idx_formPrint(Ah, row_idx, col_idx, pSpace, 0);
	}
	
	return rou_outputs;
}
