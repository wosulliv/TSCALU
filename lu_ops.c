/**
 * @file 	lu_ops.c
 * @brief 	Contains functions relating to TS matrix 
 * 		population and manipulation, as well as
 *		a suite of pivoting operations required
 *		for performing tournament pivoting.
 *		Additionally, this file contains a
 *		function for explicit LUP factorisation.
 * @author 	William O'Sullivan
 * @version 	2.0
 * @date 	2021-07-13
 */

/* IMPORTING LIBRARIES/MODULES/HEADER FILES */
#include"lu_ops.h"

/**
 * @brief      Allocates random values to an mxn matrix.
 * @params[in] matrix_memory    pointer to allocated memory for
 *                              the matrix 
 * @params[in] m_dim            dimension of the rows
 * @params[in] n_dim            dimension of the columns 
 */
int lu_build_matrix(double * const matrix_memory, const int m_dim, const int n_dim)
{
        const unsigned long seed = 170921;
        gsl_rng * r;    				/* Define generator. */
        r = gsl_rng_alloc(gsl_rng_lecuyer21); 		/* Specify type of rng; in this case L'Ecuyer. */
        gsl_rng_set(r,seed); 				/* Utilise the rng seed. */
	/* Assign random values to the matrix. */
        for(int i = 0; i < m_dim; i++){for(int j = 0; j < n_dim; j++){matrix_memory[i*n_dim+j] = gsl_ran_gaussian(r,1);}}
        return 1; 					/* return value for error checking */
}

/**
 * @brief      Prints TS matrix to terminal
 * @params[in] matrix_memory    pointer to allocated memory for
 * 				the matrix
 * @params[in] m_dim		dimension of the rows
 * @params[in] n_dim		dimension of the columns	
 */
int lu_print_matrix(double * const matrix_memory, const int m_dim, const int n_dim)
{
        fprintf(stdout, "\n \n");
        for(int i = 0; i < m_dim; i++)
        {
                fprintf(stdout, "|");
                for(int j = 0; j < n_dim; j++)
                {
                        fprintf(stdout, " %07.2f |", matrix_memory[i*n_dim + j]);
                }
                fprintf(stdout, "\n \n");
        }
        return 1; 					/* return value for err checking */
}

/**
 * @brief       For performing matrix multiplication P*A
 * @params[in]  PA      Product of matrices P and A.    
 * @params[in]  P       Matrix of size mxm.
 * @params[in]  A       Matrix of size mxn.
 * @params[in]  m_dim   Number of rows in the PA, P and A matrices.
 * @params[in]  n_dim   Number of columns in the PA and A matrices.
 */
void calculate_PA(double * const PA, double * const P, double * const A, const int m_dim, const int n_dim)
{
        double alpha = 1.0; double beta = 0.0;
        /* PA[mxn]-shaped dgemm call */ 
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m_dim, n_dim, m_dim, alpha, P, m_dim, A, n_dim, beta, PA, n_dim);
}


/**
 * @brief 	Performs explicit LUP factorisation. 
 * @params[in]  A	Matrix to be factorised.
 * @params[in]	L	Lower matrix to be populated.
 * @params[in] 	U 	Upper matrix to be populated.
 * @params[in]	P	Pivot matrix to be populated.
 * @params[in]	m_dim	Number of rows in the A, L, U and P matrices.
 * @params[in] 	n_dim	Number of columns in the A, L and U matrices.
 */
void LUP(double * const A, double * const L, double * const U, double * const P, const int m_dim, const int n_dim)
{
	/* INITIALISE MATRICES */
	for(int i = 0; i < m_dim; i++){for(int j = 0; j < n_dim; j++){U[i*n_dim+j] = A[i*n_dim+j]; if(i==j) {L[i*n_dim+j] = 1.0;}}}
	for(int i = 0; i < m_dim; i++){for(int j = 0; j < m_dim; j++){if(i==j) {P[i*m_dim+j] = 1.0;}}}

	
	/* MAIN LOOP */
	for(int k = 0; k < n_dim; k++) /* First design choice - Trefethen & Bau specified from k=1 to m-1, so in C we use k=0 to n for an mxn matrix */
	{
		/* Select i >= k to maximise U[i*n_dim+k] */
		double tmp = 0; int index = 0; for(int i = k; i < m_dim; i++){if(fabs(U[i*n_dim+k])>tmp){tmp = fabs(U[i*n_dim+k]); index = i;}}
	
		/* SWAP U ROWS */
		for(int kay = k; kay < n_dim; kay++){double z_U = U[k*n_dim+kay]; U[k*n_dim+kay] = U[index*n_dim+kay]; U[index*n_dim+kay] = z_U;}	

		/* SWAP L ROWS */
		for(int kay = 0; kay < k; kay++){double z_L = L[k*n_dim+kay]; L[k*n_dim+kay] = L[index*n_dim+kay]; L[index*n_dim+kay] = z_L;}

		/* SWAP P ROWS */
		//for(int kay = 0; kay < n_dim; kay++){double z_P = P[k*n_dim+kay]; P[k*n_dim+kay] = P[index*n_dim+kay]; P[index*n_dim+kay] = z_P;} /* Deprecated: Used to form P[mxn] */
		for(int kay = 0; kay < m_dim; kay++){double z_P = P[k*m_dim+kay]; P[k*m_dim+kay] = P[index*m_dim+kay]; P[index*m_dim+kay] = z_P;}

		/* PERFORM GAUSSIAN ELIMINATION */
		for(int j = k+1; j < m_dim; j++){L[j*n_dim+k] = U[j*n_dim+k]/U[k*n_dim+k]; 
		for(int kay = k; kay < n_dim; kay++) {U[j*n_dim+kay] -= L[j*n_dim+k] * U[k*n_dim+kay];}}
		
	}
}	

/**
 * @brief  	Performs pivoting on a sub-matrix, tracking the permutation of row indices by forming k-i pairs.
 * @params[in]	A       Sub-matrix to undergo pivoting.
 * @params[in]  m_dim   Number of rows in the A matrix.
 * @params[in]  n_dim   Number of columns in the A matrix.
 * @params[in] 	local_index	Array of indices corresponding to A.
 * @params[in]	local_k	Array for storing k values.
 * @params[in]	local_i	Array for storing i values.
 */
void pivot_contestants(double * const A, const int m_dim, const int n_dim, int * const local_index, int * const local_k, int * const local_i)
{
	for(int k = 0; k < n_dim; k++)
        {
		/* Select i >= k to maximise A[i*n_dim+k] */
		double tmp = 0; int index = 0; for(int i = k; i < m_dim; i++){if(fabs(A[i*n_dim+k])>tmp){tmp = fabs(A[i*n_dim+k]); index = i;}}
		/* k-i pairs act as the unit of pivoting, denoting the two rows being interchanged. */
		local_k[k] = local_index[k];		/* Assign the values of k. */
		local_i[k] = local_index[index];	/* Assign the values of i. */
		/* Swap the rows of the matrix to represent the permutation. */
		for(int kay = 0; kay < n_dim; kay++){double z_A = A[k*n_dim+kay]; A[k*n_dim+kay] = A[index*n_dim+kay]; A[index*n_dim+kay] = z_A;}
		/* Swap the index rows to represent and track the permutation. */
		int ind_tmp = local_index[k]; local_index[k] = local_index[index]; local_index[index] = ind_tmp;
	}
}

/**
 * @brief 	Recovers the top m_dim rows of sub-matrix A. 
 * @params[in] 	A 	Matrix of size ?xn.
 * @params[in]  m_dim   Depth of collection for a round of tournament; elsewhere, written as b.
 * @params[in]  n_dim   Number of columns in the A and W matrix.
 * @params[out]	W	Matrix of size mxn.
 */
double * collect_winners(double * const A, const int m_dim, const int n_dim)
{
	double *W = mkl_malloc(sizeof *W * (m_dim*n_dim), 64);
	for(int i = 0; i < m_dim; i++){for(int j = 0; j < n_dim; j++){W[i*n_dim+j] = A[i*n_dim+j];}}
	return W;
}

// overwrite A with new matrix!
/**
 * @brief       For performing matrix multiplication P*A
 * @params[in]  A      	Sub-matrix memory to be overwritten.    
 * @params[in]  W       Matrix of winners from the receiving processor [m/2 x n].
 * @params[in]  W_buff  Matrix of winners from the sending processor [m/2 x n].
 * @params[in]  m_dim   Number of rows in the A matrix.
 * @params[in]  n_dim   Number of columns in the A, W and W_buff matrices.
 *
 */
void combine_winners(double * const A, double * const W, double * const W_buff, const int m_dim, const int n_dim)
{
	/* The entries from W are written to the top half of A. */
	for(int i = 0; i < m_dim/2; i++)
	{
		for(int j = 0; j < n_dim; j++)
		{
			A[i*n_dim+j] = W[i*n_dim+j];
		}
	}

	/* The entries from W_buff are written to the bottom half of A. */
	for(int i = m_dim/2; i < m_dim; i++)
	{
		for(int j = 0; j < n_dim; j++)
		{
			A[i*n_dim+j] = W_buff[(i-(m_dim/2))*n_dim+j];
		}
	}
}

/**
 * @brief  	Takes the TS matrix mat, evaluates the permutations as 
 * 		described in global_k and global_i, and applies the 
 * 		permutations to produce a tournament pivoted mat[mxn].
 * @params[in]	mat		TS matrix to undergo pivoting.
 * @params[in]  global_index    Array of indices for the matrix mat.
 * @params[in] 	global_k	Array of global k indices collected from tournament pivoting.
 * @params[in] 	global_i	Array of global i indices collected from tournament pivoting.
 * @params[in] 	local_k		Array of local k indices (the final indices created on root).
 * @params[in]	local_i		Array of local i indices (the final indices created on root).
 * @params[in]  m_dim   	Number of rows in the mat matrix.
 * @params[in]  n_dim   	Number of columns in the mat and sub-matrices.
 * @params[in]  block_height	Number of rows in the sub-matrices.
 * @params[in]  log2size	Number of rounds of tournament pivoting.
 * @params[in]  size		Size of MPI communicator.
 */
double * perform_pivoting(double * const mat, int * const global_index, int * const global_k, int * const global_i, int * const local_k, int * const local_i, const int m_dim, const int n_dim, const int log2size, const int size)
{
	double *mat_tmp = mkl_malloc(sizeof *mat_tmp * (m_dim * n_dim), 64);

	int *global_index_tmp = mkl_malloc(sizeof *global_index_tmp * (m_dim), 64);
	for(int i = 0; i < m_dim; i++){global_index_tmp[i] = i;}
	
	/* global_k and global_i are designed to store all of the k-i pairs generated
 	 * at each round of the tournament. When a processor is inactive according to
 	 * the BRT, the entries from that processor are all set to -1.
 	 * As such, we have a set of arrays that store k-i pairs and preserve the 
 	 * pattern of data storage from each processor.	
 	 */

	/* The true permutations are calculated by parsing the entries of global_k and global_i to permute the global_index. */
	for(int r = 0; r < log2size*size*n_dim; r++)
	{
		/* We pass over any k-i entries (-1,-1) as this means no pivoting took place on that specific processor. */
		if((global_k[r] != -1) && (global_i[r] != -1))
		{
			if(global_k[r] == global_i[r]){continue;} /* We also pass over any k-i entries where the values are identical, i.e. there is no actual pivot. */	
			/* k-i pairs specify the IDs of rows to be permuted; swap_indices allows for rows to be permuted according to these k-i pairs. */
			swap_indices(global_index, global_index_tmp, global_k[r], global_i[r]);
		}
	}
	
	/* The true permutations also require the k-i pairs determined by the permutations of the final sub-matrix on rank 0. */
	for(int l = 0; l < n_dim; l++)
	{
		if(local_k[l] == local_i[l]){continue;} /* As before, we ignore k-i pairs where no actual pivot occurs. */
		/* We must also use the local_k and local_i created by the final round of tournament pivoting which takes place on rank 0. */
		swap_indices(global_index, global_index_tmp, local_k[l], local_i[l]);
	}


	/* We overwrite the values of mat with the values from mat_tmp as specified by the change in global_index. */
	/* NOTE: We could actually save on writes to memory by leaving mat_tmp uninitialised and assigning the values from mat to mat_tmp and returning mat_tmp instead of void. */
	for(int i = 0; i < m_dim; i++)
	{
		for(int kay = 0; kay < n_dim; kay++)
                {
                	mat_tmp[i*n_dim+kay] = mat[global_index[i]*n_dim+kay];
		}
	}

	/* Free temporary memory. */
	mkl_free(global_index_tmp);

	return mat_tmp;
}

/**
 * @brief       Takes the TS matrix mat, evaluates the permutations as 
 *              described in global_k and global_i. *DOES NOT* apply the 
 *              permutations to produce a pivoted mat[mxn].
 * @params[in]  global_index    Array of indices for the matrix mat.
 * @params[in]  global_k        Array of global k indices collected from tournament pivoting.
 * @params[in]  global_i        Array of global i indices collected from tournament pivoting.
 * @params[in]  local_k         Array of local k indices (the final indices created on root).
 * @params[in]  local_i         Array of local i indices (the final indices created on root).
 * @params[in]  m_dim           Number of rows in the mat matrix.
 * @params[in]  n_dim           Number of columns in the mat and sub-matrices.
 * @params[in]  block_height    Number of rows in the sub-matrices.
 * @params[in]  log2size        Number of rounds of tournament pivoting.
 * @params[in]  size            Size of MPI communicator.
 */
void perform_pivoting2(int * const global_index, int * const global_k, int * const global_i, int * const local_k, int * const local_i, const int m_dim, const int n_dim, const int log2size, const int size)
{
        int *global_index_tmp = mkl_malloc(sizeof *global_index_tmp * (m_dim), 64);
        for(int i = 0; i < m_dim; i++){global_index_tmp[i] = i;}

        /* global_k and global_i are designed to store all of the k-i pairs generated
         * at each round of the tournament. When a processor is inactive according to
         * the BRT, the entries from that processor are all set to -1.
         * As such, we have a set of arrays that store k-i pairs and preserve the 
         * pattern of data storage from each processor. 
         */

        /* The true permutations are calculated by parsing the entries of global_k and global_i to permute the global_index. */
        for(int r = 0; r < log2size*size*n_dim; r++)
        {
                /* We pass over any k-i entries (-1,-1) as this means no pivoting took place on that specific processor. */
                if((global_k[r] != -1) && (global_i[r] != -1))
                {
                        if(global_k[r] == global_i[r]){continue;} /* We also pass over any k-i entries where the values are identical, i.e. there is no actual pivot. */
                        /* k-i pairs specify the IDs of rows to be permuted; swap_indices allows for rows to be permuted according to these k-i pairs. */
                        swap_indices(global_index, global_index_tmp, global_k[r], global_i[r]);
                }
        }

        /* The true permutations also require the k-i pairs determined by the permutations of the final sub-matrix on rank 0. */
        for(int l = 0; l < n_dim; l++)
        {
                if(local_k[l] == local_i[l]){continue;} /* As before, we ignore k-i pairs where no actual pivot occurs. */
                /* We must also use the local_k and local_i created by the final round of tournament pivoting which takes place on rank 0. */
                swap_indices(global_index, global_index_tmp, local_k[l], local_i[l]);
        }

        /* Free temporary memory. */
        mkl_free(global_index_tmp);
}

/**
 * @brief       Performs pivoting using local indices as yielded by
 * 		result of scatter operation after invocation of
 * 		"Perform_Pivoting2"
 * @params[in]  mat			Original matrix memory [mxn]
 * @params[in]  local_index		Indices prepared for pivoting
 * @params[in]  block_row_memory	Memory of panel sub-matrix
 * @params[in]  m_dim			m dimension of panel
 * @params[in]  n_dim			n dimension of panel
 */
void permute_matrix(double * const mat, int * const local_index, double * const block_row_memory, const int m_dim, const int n_dim)
{
        for(int i = 0; i < m_dim; i++)
        {
                for(int kay = 0; kay < n_dim; kay++)
                {
                        block_row_memory[i*n_dim+kay] = mat[local_index[i]*n_dim+kay];
                }
        }
}


/**
 * @brief  	Keeps track of and performs the row swapping operations.
 * @params[in] 	global_index		Global index which is permuted according to the given k-i pair.
 * @params[in]	global_index_tmp	Temporary global index which records the locations of rows according to k-i pairs.
 * @params[in]	k			k ID of row to be swapped.
 * @params[in] 	i 			i ID of row to be swapped.
 */
void swap_indices(int * const global_index, int * const global_index_tmp, int k, int i)
{
	/* Take the indices of permutation from global_index_tmp*/
	int k_tmp = global_index_tmp[k];
	int i_tmp = global_index_tmp[i];

	/* Swap the indices stored in the global index according to the values pulled from global_index_tmp. */
	int ind_tmp = global_index[k_tmp];
        global_index[k_tmp] = global_index[i_tmp];
        global_index[i_tmp] = ind_tmp;	

	/* Swap the indices in global_index_tmp to reflect the pivoting of k and i. */
	global_index_tmp[k] = i_tmp;
	global_index_tmp[i] = k_tmp;
}

/**
 * @brief       Performs gaussian elimination to form L and U matrices. 
 * @params[in]  A       Matrix to undergo elimination.
 * @params[in]  L       Lower matrix to be populated.
 * @params[in]  U       Upper matrix to be populated.
 * @params[in]  m_dim   Number of rows in the A, L, and U matrices.
 * @params[in]  n_dim   Number of columns in the A, L and U matrices.
 */
void gaussian_elimination(double * const A, double * const L, double * const U, const int m_dim, const int n_dim)
{
        /* Performs initialisation of L and U */
        for(int i = 0; i < m_dim; i++){for(int j = 0; j < n_dim; j++){U[i*n_dim+j] = A[i*n_dim+j]; if(i==j) {L[i*n_dim+j] = 1.0;}}}
        /* Loops through columns of L and U to perform Gaussian Elimination. */
        for(int k = 0; k < n_dim; k++){
        for(int j = k+1; j < m_dim; j++){L[j*n_dim+k] = U[j*n_dim+k]/U[k*n_dim+k];
        for(int kay = k; kay < n_dim; kay++) {U[j*n_dim+kay] -= L[j*n_dim+k] * U[k*n_dim+kay];}}}
}

/**
 * @brief       Reshapes U matrix for multiplication [nxn].
 * @params[in]  U       Upper matrix to be reshaped.
 * @params[in]  m_dim   Number of rows in the U matrix.
 * @params[in]  n_dim   Number of columns in the U matrix, and rows and columns of U_new.
 */
double * reshape_U(double * const U, const int n_dim)
{
	double * U_new = mkl_malloc(sizeof(*U_new)*(n_dim*n_dim), 64);
	for(int i = 0; i < n_dim; i++){for(int j = 0; j < n_dim; j++){U_new[i*n_dim+j] = U[i*n_dim+j];}}
	//lu_print_matrix(U_new, n_dim, m_dim);
	return U_new;
}

/**
 * @brief       Reshapes L matrix for multiplication [mxm].
 * @params[in]  L       Lower matrix to be reshaped.
 * @params[in]  m_dim   Number of rows in the L matrix, and the rows and columns of L_new.
 * @params[in]  n_dim   Number of columns in the L matrix.
 */
double * reshape_L(double * const L, const int m_dim, const int n_dim)
{
        double * L_new = mkl_malloc(sizeof(*L_new)*(m_dim*m_dim), 64);
	for(int i = 0; i < m_dim; i++){for(int j = 0; j < m_dim; j++){if(i==j){L_new[i*m_dim+j] = 1.0;}else{L_new[i*m_dim+j] = 0.0;}}}
	for(int i = 0; i < m_dim; i++)
	{
		for(int j = 0; j < n_dim; j++)
		{
			L_new[i*(m_dim/n_dim)*n_dim+j] = L[i*n_dim+j];
		}
	}
        //lu_print_matrix(L_new, m_dim, m_dim);
        return L_new;
}

/**
 * @brief       For performing matrix multiplication L*U
 * @params[in]  LU      Product of matrices L and U.    
 * @params[in]  L       Matrix of size mxm.
 * @params[in]  U       Matrix of size nxn.
 * @params[in]  m_dim   Number of rows in the LU and L matrices.
 * @params[in]  n_dim   Number of columns in the PA and U matrices.
 */
void calculate_LU(double * const LU, double * const L, double * const U, const int m_dim, const int n_dim)
{
        double alpha = 1.0; double beta = 0.0;
        /* LU[mxn]-shaped dgemm call */ 
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, m_dim, n_dim, m_dim, alpha, L, m_dim, U, n_dim, beta, LU, n_dim);
}


/**
 * @brief       Calculates the infinity norm of a matrix.
 * @params[in]  matrix  Target matrix of inf norm calculation.    
 * @params[in]  m_dim   Number of rows in the matrix.
 * @params[in]  n_dim   Number of columns in the matrix.
 */
double infinity_norm(double * const matrix, const int m_dim, const int n_dim)
{
	double max = 0;
	for(int i = 0; i < m_dim; i++)
	{
		double sum = 0;
		for(int j = 0; j < n_dim; j++)
		{
			sum += fabs(matrix[i*n_dim+j]);
		}
		if(sum > max)
		{
			max=sum;
		}
	}
	return max;
}


