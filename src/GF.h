#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>
#include <vector>
#include <cmath>
using namespace std;

// Print 1D ptr array based in the default order 
template <typename TypeA>
void GF_rawPrint(TypeA *ptr1D, int size)
{
	for (int i = 0; i < size; i++)
		cout << ptr1D[i] << " ";
	cout << endl << endl;	
}

// Print 2D ptr array based on the default row/column order 
template <typename TypeA>
void GF_rawPrint(TypeA **ptr2D, int rows, int cols)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
			cout << ptr2D[i][j] << " ";
		cout << endl;
	}
	cout << endl;
}

// Print 1D ptr array based in the specified order 
template <typename TypeA>
void GF_idx_rawPrint(TypeA *ptr1D, vector<int> &idx)
{
	for (int i = 0; i < idx.size(); i++)
		cout << ptr1D[idx[i]] << "\n";
	cout << endl << endl;	
}

// Print 2D ptr array based on a specific row/column order 
template <typename TypeA>
void GF_idx_rawPrint(TypeA **ptr2D, vector<int> &row_idx, vector<int> &col_idx)
{
	for (int i = 0; i < row_idx.size(); i++)
	{
		for (int j = 0; j < col_idx.size(); j++)
			cout << ptr2D[row_idx[i]][col_idx[j]] << " ";
		cout << endl;
	}
	cout << endl;
}

// Print 2D mpz_t ptr array based on a specific row/column order 
void GFz_idx_formPrint(mpz_t **ptr2D, vector<int> &row_idx, vector<int> &col_idx, int pSpace, int col_adj)
{
	double temp;
	int stringSize_entry;

	int pSpace_adj[col_idx.size()];

	for (int j = 0; j < col_idx.size(); j++)
	{
		pSpace_adj[j] = min(pSpace + col_adj, pSpace + j);
		if(j >= 10 && col_adj>= 8)
			pSpace_adj[j] += 2;
	}
	for (int i = 0; i < row_idx.size(); i++)
	{
		for (int j = 0; j < col_idx.size(); j++)
		{			
			temp =  mpz_get_d( ptr2D[row_idx[i]][col_idx[j]]);
			
			if(temp == 0)
				stringSize_entry = 1;
			else
				stringSize_entry = ceil(log10(std::abs((double)temp)+1));

			if(temp >= 0)
				cout << " ";

			cout << ptr2D[row_idx[i]][col_idx[j]];

			for (int k=0; k < pSpace_adj[j] - stringSize_entry; k++)
				cout << " ";

			cout << " ";
		}
		cout << endl;
	}
	cout << endl;
}

// Print 2D mpz_t ptr array based on a specific row/column order 
void GFz_idx_framePrint(mpz_t **ptr2D, vector<int> &row_idx, vector<int> &col_idx, vector<int> &frame_sgn, int pSpace, int col_adj)
{
	mpz_t posEntry;
	mpz_t mpz_ZERO;
	mpz_inits(mpz_ZERO, posEntry, NULL);

	double temp;
	int stringSize_entry;

	int pSpace_adj[col_idx.size()];

	for (int j = 0; j < col_idx.size(); j++)
	{
		pSpace_adj[j] = min(pSpace + col_adj, pSpace + j);
		if (j >= 10 && col_adj >= 8)
			pSpace_adj[j] += 2;
	}

	for (int i = 0; i < row_idx.size(); i++)
	{
		//Frame index is related to the column index
		for (int j = 0; j < i; j++)
		{
			temp = mpz_get_d(ptr2D[row_idx[i]][col_idx[j]]);

			if (temp == 0)
				stringSize_entry = 1;
			else
				stringSize_entry = ceil(log10(std::abs((double)temp) + 1));

			if (temp >= 0)
				cout << " ";

			//Print current entries as is
			if (frame_sgn[j] == 1)
				cout << ptr2D[row_idx[i]][col_idx[j]];

			else
			{
				if (mpz_cmp(ptr2D[row_idx[i]][col_idx[j]], mpz_ZERO) <= 0)
				{
					mpz_abs(posEntry, ptr2D[row_idx[i]][col_idx[j]]);
					cout << posEntry;
				}

				else
				{
					cout << "-" << ptr2D[row_idx[i]][col_idx[j]];
				}
			}

			for (int k = 0; k < pSpace_adj[j] - stringSize_entry; k++)
				cout << " ";

			cout << " ";
		}


		//Frame index is related to the row index
		//Print current entries as is
		if (frame_sgn[i] == 1)
		{
			for (int j = i; j < col_idx.size(); j++)
			{
				temp = mpz_get_d(ptr2D[row_idx[i]][col_idx[j]]);

				if (temp == 0)
					stringSize_entry = 1;
				else
					stringSize_entry = ceil(log10(std::abs((double)temp) + 1));

				if (temp >= 0)
					cout << " ";

				cout << ptr2D[row_idx[i]][col_idx[j]];

				for (int k = 0; k < pSpace_adj[j] - stringSize_entry; k++)
					cout << " ";

				cout << " ";
			}
		}

		//Print opposite sign for each entry
		else
		{
			for (int j = i; j < col_idx.size(); j++)
			{
				temp = mpz_get_d(ptr2D[row_idx[i]][col_idx[j]]);

				if (temp == 0)
					stringSize_entry = 1;
				else
					stringSize_entry = ceil(log10(std::abs((double)temp) + 1));

				if (temp >= 0)
					cout << " ";

				if (mpz_cmp(ptr2D[row_idx[i]][col_idx[j]], mpz_ZERO) <= 0)
				{
					mpz_abs(posEntry, ptr2D[row_idx[i]][col_idx[j]]);
					cout << posEntry;
				}

				else
				{
					cout << "-" << ptr2D[row_idx[i]][col_idx[j]];
				}
				for (int k = 0; k < pSpace_adj[j] - stringSize_entry; k++)
					cout << " ";

				cout << " ";
			}
		}
		cout << endl;
	}
	cout << endl;
}

// Print 2D mpq_t ptr array based on a specific row/column order 
void GFq_idx_formPrint(mpq_t **ptr2D, vector<int> &row_idx, vector<int> &col_idx, int pSpace, int col_adj)
{
	double temp;
	int stringSize_entry;

	int pSpace_adj[col_idx.size()];

	for (int j = 0; j < col_idx.size(); j++)
	{
		pSpace_adj[j] = min(pSpace + col_adj, pSpace + j);
		if(j >= 10 && col_adj>= 8)
			pSpace_adj[j] += 2;
	}
	for (int i = 0; i < row_idx.size(); i++)
	{
		for (int j = 0; j < col_idx.size(); j++)
		{			
			temp =  mpq_get_d( ptr2D[row_idx[i]][col_idx[j]]);
			
			if(temp == 0)
				stringSize_entry = 1;
			else
				stringSize_entry = ceil(log10(abs(temp)+1));

			if(temp >= 0)
				cout << " ";

			cout << ptr2D[row_idx[i]][col_idx[j]];

			for (int k=0; k < pSpace_adj[j] - stringSize_entry; k++)
				cout << " ";

			cout << " ";
		}
		cout << endl;
	}
	cout << endl;
}

// Initalize an mpz matrix of the specified dimensions to 0
mpz_t** GFz_mat_init_zeros(int rows, int cols)
{
	mpz_t ** mpz_mat;
	mpz_mat = new mpz_t *[rows];			
	for (int i = 0; i < rows; i++)
		mpz_mat[i] = new mpz_t[cols];

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
			mpz_init(mpz_mat[i][j]);
	}
	return mpz_mat;
}

// Initalize an mpq matrix of the specified dimensions to 0
mpq_t** GFq_mat_init_zeros(int rows, int cols)
{
	mpq_t ** mpq_mat;
	mpq_mat = new mpq_t *[rows];			
	for (int i = 0; i < rows; i++)
		mpq_mat[i] = new mpq_t[cols];

	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
			mpq_init(mpq_mat[i][j]);
	}
	return mpq_mat;
}

// Initalize an mpz vector of the specified dimension to 0
mpz_t* GFz_vec_init_zeros(int rows)
{
	mpz_t * mpz_vec;
	mpz_vec = new mpz_t[rows];
	for (int i = 0; i < rows; i++)
		mpz_init(mpz_vec[i]);
	return mpz_vec;
}

// Initalize an mpq vector of the specified dimension to 0
mpq_t* GFq_vec_init_zeros(int rows)
{
	mpq_t * mpq_vec;
	mpq_vec = new mpq_t[rows];
	for (int i = 0; i < rows; i++)
		mpq_init(mpq_vec[i]);
	return mpq_vec;
}

// Initalize a randomized mpz matrix according to the specified arguments
mpz_t** GFz_mat_init_density(int rows, int cols, int lb, int ub, int seed1, int seed2, double density, bool diag)
{	
	mpz_t ** mpz_mat = GFz_mat_init_zeros(rows, cols);

	mpz_t temp_z, mpz_ZERO;
	mpz_inits(temp_z, mpz_ZERO, NULL);

	int idx1, idx2;
	int rn; 
	double fraction_positive;

	int numNonZeros = ceil(density*rows*cols);

	//To determine the right proportion of positive to negative entries in input matrix
	if (lb >= 0)
		fraction_positive = 1;
	else if(ub <= 0)
		fraction_positive = 0;
	else
		fraction_positive = (double) ub/(ub-lb);


	//Make sure the (rows)x(rows) initial left-hand subArix is nonsingular
	//This is done by placing a nonzero element along each column without 
	//repeating rows; do this so the location is randomized. 

	srand(seed1);
	vector<int> row_idx_rand;

	for(int i=0; i < rows; i++)
		row_idx_rand.push_back(i);

	//This creates an initial non-singular left-hand Arix that likely has diagonal 
	//elements equal to zero
	if(!diag)
	{
		for(int i=0; i < 2*rows; i++)
		{
			idx1 = rand() % rows;
			idx2 = rand() % rows;
			swap(row_idx_rand[idx1],row_idx_rand[idx2]);
		}
	}

	srand(2*seed1);

	int temp_count = 0;

	//Mark placement of initial nonzero values to guarantee nonsingularity when mpz_mat is 
	//a square matrix (or has more columns than rows)
	for(int i=0; i < min(rows,cols); i++)
	{
		idx1 = row_idx_rand[i];
		if((double) rand() / (RAND_MAX) < fraction_positive)
			mpz_set_si(mpz_mat[idx1][i], 1);
		else
			mpz_set_si(mpz_mat[idx1][i], -1);
		temp_count++;
	}

	numNonZeros -= temp_count;

	//Mark placement of remaining initial nonzero values
	while(numNonZeros > 0)
	{
		idx1 = rand() % rows;
		idx2 = rand() % cols;

		if(mpz_cmp(mpz_mat[idx1][idx2],mpz_ZERO)==0)
		{
			if((double) rand() / (RAND_MAX) < fraction_positive)
				mpz_set_si(mpz_mat[idx1][idx2], 1);
			else
				mpz_set_si(mpz_mat[idx1][idx2], -1);
			numNonZeros--;
		}
	}

	//For random number generation
	srand(seed2);

	//Substitute nonzero tokens with random nonzero values 
	for (int j = 0; j < cols; j++)
	{
		for (int i = 0; i < rows; i++)
		{
			if(mpz_cmp(mpz_mat[i][j],mpz_ZERO)<0)
			{
				mpz_set_si(mpz_mat[i][j], -1 * (rand() % (-1*lb) + 1));
			}
			else if(mpz_cmp(mpz_mat[i][j],mpz_ZERO)>0)
				mpz_set_si(mpz_mat[i][j], rand() % ub + 1);
		}
	}	
	return mpz_mat;
}

// Initalize a randomized mpq matrix according to the specified arguments
mpq_t** GFq_mat_init_density(int rows, int cols, int lb, int ub, int seed1, int seed2, double density, bool diag)
{	
	mpq_t ** mpq_mat = GFq_mat_init_zeros(rows, cols);

	mpq_t temp_z, mpq_ZERO;
	mpq_inits(temp_z, mpq_ZERO, NULL);

	int idx1, idx2;
	int rn; 
	double fraction_positive;

	int numNonZeros = ceil(density*rows*cols);

	//To determine the right proportion of positive to negative entries in input Arix
	if (lb >= 0)
		fraction_positive = 1;
	else if(ub <= 0)
		fraction_positive = 0;
	else
		fraction_positive = (double) ub/(ub-lb);


	//Make sure the (rows)x(rows) initial left-hand subArix is nonsingular
	//This is done by placing a nonzero element along each column without 
	//repeating rows; do this so the location is randomized. 

	srand(seed1);
	vector<int> row_idx_rand;

	for(int i=0; i < rows; i++)
		row_idx_rand.push_back(i);

	//This creates an initial non-singular left-hand Arix that likely has diagonal 
	//elements equal to zero
	if(!diag)
	{
		for(int i=0; i < 2*rows; i++)
		{
			idx1 = rand() % rows;
			idx2 = rand() % rows;
			swap(row_idx_rand[idx1],row_idx_rand[idx2]);
		}
	}

	srand(2*seed1);

	int temp_count = 0;

	//Mark placement of initial nonzero values to guarantee nonsingularity when mpq_mat is 
	//a square matrix (or has more columns than rows)
	for(int i=0; i < min(rows,cols); i++)
	{
		idx1 = row_idx_rand[i];
		if((double) rand() / (RAND_MAX) < fraction_positive)
			mpq_set_si(mpq_mat[idx1][i], 1, 1);
		else
			mpq_set_si(mpq_mat[idx1][i], -1, 1);
		temp_count++;
	}

	numNonZeros -= temp_count;

	//Mark placement of remaining initial nonzero values
	while(numNonZeros > 0)
	{
		idx1 = rand() % rows;
		idx2 = rand() % cols;

		if(mpq_cmp(mpq_mat[idx1][idx2],mpq_ZERO)==0)
		{
			if((double) rand() / (RAND_MAX) < fraction_positive)
				mpq_set_si(mpq_mat[idx1][idx2], 1, 1);
			else
				mpq_set_si(mpq_mat[idx1][idx2], -1, 1);
			numNonZeros--;
		}
	}

	//For random number generation
	srand(seed2);

	//Substitute nonzero tokens with random nonzero values 
	for (int j = 0; j < cols; j++)
	{
		for (int i = 0; i < rows; i++)
		{
			if(mpq_cmp(mpq_mat[i][j],mpq_ZERO)<0)
			{
				mpq_set_si(mpq_mat[i][j], -1 * (rand() % (-1*lb) + 1), 1);
			}
			else if(mpq_cmp(mpq_mat[i][j],mpq_ZERO)>0)
				mpq_set_si(mpq_mat[i][j], rand() % ub + 1, 1);
		}
	}	
	return mpq_mat;
}

// Initalize a randomized mpz vector according to the specified arguments
mpz_t* GFz_vec_init_rand(int rows, int lb, int ub, int seed)
{
	mpz_t * mpz_vec = GFz_vec_init_zeros(rows);

	//For random number generation
	srand(seed);

	for (int i = 0; i < rows; i++)
		mpz_set_si(mpz_vec[i], rand() % ub + 1);

	return mpz_vec;
}

// For copying mpq matrices
void GFq_mat_copy(mpq_t **mat_orig, mpq_t **mat_copy, int rows, int cols)
{
	for (int i=0; i < rows; i++)
	{
		for (int j=0; j < cols; j++)
		{
			mpq_init(mat_copy[i][j]);
			mpq_set(mat_copy[i][j], mat_orig[i][j]);
		}
	}
}

// For copying mpz matrices
void GFz_mat_copy(mpz_t **mat_orig, mpz_t **mat_copy, int rows, int cols)
{
	for (int i=0; i < rows; i++)
	{
		for (int j=0; j < cols; j++)
		{
			mpz_init(mat_copy[i][j]);
			mpz_set(mat_copy[i][j], mat_orig[i][j]);
		}
	}
}

// For copying mpz vectors 
void GFz_vec_copy(mpz_t *vec_orig, mpz_t *vec_copy, int rows)
{
	for (int i=0; i < rows; i++)
	{
		mpz_init(vec_copy[i]);
		mpz_set(vec_copy[i], vec_orig[i]);
	}
}

// For copying mpz matrices into mpz matrices
void GFzToq_mat_copy(mpz_t **mpz_mat, mpq_t **mpq_mat, int rows, int cols)
{
	for (int i=0; i < rows; i++)
	{
		for (int j=0; j < cols; j++)
		{
			mpq_init(mpq_mat[i][j]);
			mpq_set_z(mpq_mat[i][j], mpz_mat[i][j]);
		}
	}
}

// Replace a set of consecutive columns in an existing matrix with an input block of columns
void GFz_mat_insert_cols(int idx_ins, mpz_t **mpz_mat, mpz_t **mpz_col_block, vector<int> &col_idx, int rows, int cols_ins)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = idx_ins; j < idx_ins + cols_ins; j++)
			mpz_set(mpz_mat[i][col_idx[j]], mpz_col_block[i][j - idx_ins]);
	}
}

// Replace a set of consecutive columns in an existing matrix with an input block of columns
void GFq_mat_insert_cols(int idx_ins, mpq_t **mpq_mat, mpq_t **mpq_col_block, vector<int> &col_idx, int rows, int cols_ins)
{
	for (int i = 0; i < rows; i++)
	{
		for (int j = idx_ins; j < idx_ins + cols_ins; j++)
			mpq_set(mpq_mat[i][col_idx[j]], mpq_col_block[i][j-idx_ins]);
	}
}

// For printing vectors of different types
template <typename TypeA>
void GF_vecPrint(vector<TypeA> &vec)
{
	for (int i = 0; i < vec.size(); i++)
		cout << vec[i] << " ";
	cout << endl;	
}

// For setting vectors of different types to 0
template <typename TypeA>
void GF_setVecToZero(vector<TypeA> &vec, int size)
{
	vec.clear();
	for(int i = 0; i < size; i++)
		vec.push_back(0);
}

//Math Functions (self-explanatory based on the name)
template <typename TypeA>
TypeA GF_vecSum(vector<TypeA> &vec)
{
	TypeA sum = (TypeA) 0;

	for (int i = 0; i < vec.size(); i++)
		sum += vec[i];
	return sum;
}

void GFq_matmatMul(mpq_t **mat_prod, mpq_t **matL, mpq_t **matR, int rows_L, int cols_L, int cols_R)
{
	mpq_t temp;
	mpq_init(temp);
	
	mpz_t mpz_ZERO;
	mpz_init(mpz_ZERO);

	for(int i = 0; i < rows_L; i++)
	{
	    for(int j = 0; j < cols_R; j++)
	    {	
			mpq_set_z(mat_prod[i][j], mpz_ZERO);
			
	        for(int k = 0; k < cols_L; k++)
			{
				mpq_mul(temp, matL[i][k], matR[k][j]);	
				mpq_add(mat_prod[i][j], mat_prod[i][j], temp);
			}
	    }      
	}
}

void GFq_idx_matmatMul(mpq_t **mat_prod, mpq_t **matL, mpq_t **matR,  vector<int> &row_idx_L, vector<int> &col_idx_L, int cols_R)
{
	mpq_t temp;
	mpq_init(temp);
	
	mpz_t mpz_ZERO;
	mpz_init(mpz_ZERO);

	for(int i = 0; i < row_idx_L.size(); i++)
	{
	    for(int j = 0; j < cols_R; j++)
	    {	
			mpq_set_z(mat_prod[row_idx_L[i]][j], mpz_ZERO);
			
	        for(int k = 0; k < col_idx_L.size(); k++)
			{
				mpq_mul(temp, matL[row_idx_L[i]][col_idx_L[k]], matR[row_idx_L[k]][j]);	
				mpq_add(mat_prod[row_idx_L[i]][j], mat_prod[row_idx_L[i]][j], temp);
			}
	    }      
	}
}

bool GFq_matComp(mpq_t **matL, mpq_t **matR, int rows, int cols)
{
	int sum = 0;
	mpq_t temp;
	mpq_init(temp);
	
	for (int i = 0; i < rows; i++)
	{
		for (int j = 0; j < cols; j++)
			sum += mpq_equal(matL[i][j], matR[i][j]);
	}

	if(sum == rows*cols)
		return true;
	else
		return false;
}

bool GFz_vecComp(mpz_t *vec1, mpz_t *vec2, int rows)
{
	bool vecs_equal = true;
	
	for (int i = 0; i < rows; i++)
	{
		if(mpz_cmp(vec1[i],vec2[i]) != 0)
			vecs_equal = false;
	}
	return vecs_equal;
}

// Change the sign of a vector's elements starting with index k
void GF_vec_consencutive_sgn_flip(vector<int> &vec, int size, int k)
{
	for (int i = k; i < size; i++)
		vec[i] = -1 * vec[i];
}
