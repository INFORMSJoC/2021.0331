class Exact_LU_Alg
{
	int lb;
	int ub;
	int seed1;	//For determining nonzero elements of input matrix randomly
	int seed2;	//For filling in the input matrix values
	double density;

	public:

		int rows;
		int cols;
		int cols_B;
		int ex_idx;
		int up_idx_1;
		int up_idx_2;
		int num_ent;
		int num_perms;
		int stepNum;
		int pSpace; //for printing matrices in a visually appealing format

		// For mpz comparisons
		mpz_t mpz_ZERO;
		mpq_t mpq_ZERO;

		// Iterative input matrix and fixed input matrix
		mpz_t **A, **A0;
		mpq_t **Aq, **Aq0; //Rational version 

		// Matrix that will hold the rank-one updated factorization
		mpz_t **Ah;

		// Input right-hand vector(s) 
		mpz_t **B;
		mpq_t **Bq; //Rational version

		// Entering column(s)
		mpz_t **A_plus, **A_plus0;
		mpq_t **Aq_plus; //Rational version

		// Forward Substitution vector(s) 
		mpz_t **Y;
		mpq_t **Yq; //Rational version

		// Forward substitution iterative vector used in rank-one update
		mpz_t *y_k_old, *y_k, *z_k_old, *z_k;
		
		// Update vectors used in rank-one update
		mpz_t *v, *w;

		// Backward Substitution vector(s) 
		mpz_t **X_p;//This is the x'-vector (i.e., x'=det(A)x)
		mpq_t **Xq; //Rational version

		// Holds Exact solutions up to double precision
		double **Xd;

		// Standard Row and Column Index Order
		vector<int> row_idx_std;
		vector<int> col_idx_std;

		// Current Row and Column index Order
		vector<int> row_idx;
		vector<int> col_idx;


		// Current Row and Column index Order for update algorithms
		vector<int> up_row_idx;
		vector<int> up_col_idx;

		// Same but for RHS columns
		vector<int> col_idx_B;

	/*-----------------------------------------------------------------------------------------
		                        FUNCTIONS
	-----------------------------------------------------------------------------------------*/
		// Getter functions
		int get_ex_idx(void) { return ex_idx;}
		int get_num_ent(void) { return num_ent;}
		int get_lb (void){return lb;}
		int get_ub (void){return ub;}
		int get_seed1 (void){return seed1;}
		int get_seed2 (void){return seed2;}
		double get_density (void){return density;}

		// Setter functions
		void set_ex_idx(int ex_idx_) {ex_idx=ex_idx_;}
		void set_num_ent(int num_ent_) { num_ent = num_ent_;}
		void set_lb (int lb_){lb=lb_;}
		void set_ub (int ub_){ub=ub_;}
		void set_seed1 (int seed1_){seed1=seed1_;}
		void set_seed2 (int seed2_){seed2=seed2_;}
		void set_density (double density_){density=density_;}
		
		// Initializes members from cmdOpt object
		void load_from_cmdOpt(cmdOpt &run);

		// Verifies solution is correct by multiplying A0 times Xq and comparing produce with Bq 
		bool verify_sol(string);
		
		// Verifies solution is correct by comparing Bq with an alternative Bq matrix known to be exact
		bool verify_sol(string, mpq_t **);

		// Bartels Golub permutation (Modifies row and column index arrays)
		void BG_Perm();
};

//Initializing from cmdOpt object
void Exact_LU_Alg::load_from_cmdOpt(cmdOpt &run)
{
	rows = run.dim;
	cols = run.dim;
	cols_B = 0; 
	ex_idx = run.up_idx_1;
	num_ent = 0;
	lb = run.lb;
	ub = run.ub;
	seed1 = run.seed1;
	seed2 = run.seed2;
	density = run.density;

	stepNum = 0;
	pSpace = ceil(log10(get_ub()+1))+1;

	mpz_init(mpz_ZERO);
	mpq_init(mpq_ZERO);

	for (int i=0; i < rows; i++)
	{
		row_idx.push_back(i);
		row_idx_std.push_back(i);
	}

	for (int j=0; j < cols; j++)
	{
		col_idx.push_back(j);
		col_idx_std.push_back(j);
	}

	//This creates a copy of the initial matrix in rational form (to be used to verify the solution is exact)
	if (run.check_sol == 1)
		Aq0 = GFq_mat_init_density(rows, cols, get_lb(), get_ub(), get_seed1(), get_seed2(), get_density(), true);
	
}

//Verifies solution is correct by multiplying A0 times Xq and comparing product with Bq 
bool Exact_LU_Alg::verify_sol(string alg_name)
{
	mpq_t mpq_ZERO;
	mpq_init(mpq_ZERO);
	
	//Holds the matrix product: Initialize to 0
	mpq_t ** matProd =  GFq_mat_init_zeros(rows, cols_B);
			
	GFq_idx_matmatMul(matProd, Aq0, Xq, row_idx, col_idx, cols_B);

	//Cumulative difference between Aq0(Xq) and Bq
	mpq_t cumu_diff, temp;
	mpq_inits(cumu_diff, temp, NULL);

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols_B; j++)
		{
			mpq_sub(temp, matProd[row_idx[i]][j], Bq[row_idx[i]][j]);
			mpq_abs(temp, temp);
			mpq_add(cumu_diff, cumu_diff, temp);
		}
	}

	if(mpq_cmp(cumu_diff, mpq_ZERO)==0)					
		return true;
	else 
		return false;
}

//Verifies solution is correct by comparing Bq with an alternative Bq matrix known to be exact
bool Exact_LU_Alg::verify_sol(string alg_name, mpq_t **exact_Xq)
{
	return GFq_matComp(Xq, exact_Xq, rows, cols_B);
}

// Bartels Golub permutation (Modifies row and column index arrays)
void Exact_LU_Alg::BG_Perm()
{
	up_row_idx.clear();
	up_col_idx.clear();
	
	//Initialize update row and column index arrays to those used in the factorization construction
	for (int i = 0; i < rows; i++)
		up_row_idx.push_back(row_idx[i]);

	for (int j = 0; j < cols; j++)
		up_col_idx.push_back(col_idx[j]);

	for (int j = ex_idx + num_ent - 1; j >= ex_idx; j--)
	{
		for (int k = j; k < cols - 1; k++)
			std::swap(up_col_idx[k], up_col_idx[k + 1]);
	}
}
