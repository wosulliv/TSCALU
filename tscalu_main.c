/**
 * @file 	tscalu_main.c
 * @brief 	Main function for TS CALU which invokes 
 * 		communication-avoiding tournament pivoting.
 * @author 	William O'Sullivan
 * @version 	1.9
 * @date 	2021-07-19
 */

/* IMPORTING LIBRARIES/MODULES/HEADER FILES */
#include <mkl_types.h>
#include <mkl_cblas.h>
#include "mkl.h"
#include <mpi.h>
#include <unistd.h>
#include"lu_ops.h"
#include"file_ops.h"

/**
 * @brief Main function for scalable_tslu.c
 */
int main(int argc, char *argv[])
{
	/* Initialising MPI */
	MPI_Init(&argc, &argv);
	/* Beginning timer */
        MPI_Barrier(MPI_COMM_WORLD);
        double start = MPI_Wtime();

	/* Setting up MPI parameters */
	int rank, size;
        MPI_Status status;
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);

	int opt; int rows; int cols; int type;

	while((opt = getopt(argc, argv, "r:c:t:")) != -1){
		switch(opt){
			case 'r':
				rows = atoi(optarg);
				break;
			case 'c':
				cols = atoi(optarg);
				break;
			case 't':
				type = atoi(optarg);
				break;
		}
        }

	const int tag1 = 0; const int tag2 = 1; const int tag3 = 2; const int tag4 = 3;
	const int matrix_build_type = type;	/* If 0: builds a matrix [mxn] of randomly generated numbers along a normal distribution between +/- 1.
					   	 * Otherwise: takes an [mxn] panel from the matrix specified by the "filename" variable. */
	
	/* 
 	 * IMPORTANT: IF MATRIX BUILD TYPE = 0 NO MATRIX FILE WILL BE READ IN! 
 	 * ADDITIONALLY: (Only) One filename must be specified at all times! 
 	 */
	
	//const char filename[] = "matrices_4096/circulant.txt";
	//const char filename[] = "matrices_4096/diag_dom.txt";
	//const char filename[] = "matrices_4096/jordan_block.txt";
	const char filename[] = "matrices_4096/kms.txt";
	//const char filename[] = "matrices_4096/neumann.txt";

	/* Establishing matrix and sub-matrix dimensions */
	const int m_dim = rows;	/* The number of rows */
        const int n_dim = cols;	/* The number of columns */
	
	const int block_size = (m_dim / size) * n_dim;	/* The number of entries in a sub-matrix */
	const int block_height = (m_dim / size);	/* The number of rows in a sub-matrix */
	const int b = n_dim;				//block_height/2;
				 			/* The number of rows that would be passed along from each round of tournament pivoting.
							 * This value was chosen to be block_height/2 in order to simplify the use of memory at every
							 * round of the tournament. */
	const int tournament_height = 2*b;
	
	const int log2size = log2(size);              	/* This is used to specify the number of rounds of pivoting which take place */
        const int lastpower = 1 << log2size;    	/* This operation was used when the Binary Reduction Tree was being developed for a non-(2^n) number of processors.
						 	 * It has since been deprecated, but is being left in the source as a ground for future expansion.
						 	 * Bitshift Operation = 1*(2^log2size); */
	
	/* This is some preliminary bounds-testing to help ensure that the matrix given satisfies some basic requirements */
	if(rank==0){
		if(block_height%2!=0 && m_dim!=size)
		{fprintf(stderr, "\n==> ERROR: nprocs used fails to divide evenly into the column length.\n"); MPI_Abort(MPI_COMM_WORLD, 1);}
		if(size>m_dim)
		{fprintf(stderr, "\n==> ERROR: nprocs used exceeds column length.\n"); MPI_Abort(MPI_COMM_WORLD, 1);}
		if((m_dim/n_dim)<size)
		{fprintf(stderr, "\n==> ERROR: The matrix is being decomposed into matrices of the wrong dimenstion; n*m rather than m*n.\n");MPI_Abort(MPI_COMM_WORLD, 1);}
		if(((m_dim/size)/2)<n_dim)
                {fprintf(stderr, "\n==> ERROR: The matrix is too wide to be decomposed on this core count.\n");MPI_Abort(MPI_COMM_WORLD, 1);}}	

	double start_malloc = MPI_Wtime();
	/* Allocating memory for and building the initial matrix, as well as allocating memory for the L and U matrices. */
	double *mat = mkl_malloc(sizeof *mat * (m_dim * n_dim), 64);
	if(matrix_build_type == 0){lu_build_matrix(mat,m_dim,n_dim);}
	else{read_panel(mat, m_dim, n_dim, filename);}
	double *A = mkl_malloc(sizeof *A * (block_height*n_dim), 64);
	for(int i = 0; i < block_height; i++){for(int j = 0; j < n_dim; j++){A[i*n_dim+j] = mat[(i*n_dim+j) + rank*block_size];}}
	double *A_minor = mkl_malloc(sizeof *A_minor * (tournament_height*n_dim), 64);
	double *L = mkl_malloc(sizeof *L * (m_dim*n_dim), 64);
	double *U = mkl_malloc(sizeof *U * (m_dim*n_dim), 64);
	double *LU = mkl_malloc(sizeof *LU * (m_dim*n_dim), 64);

	double original_mat_inf_norm = infinity_norm(mat, m_dim, n_dim);

	/* Allocating memory for and building index systems - one for global entries, and one for local entries. */
	int *global_index = mkl_malloc(sizeof *global_index * (m_dim), 64);
	for(int i = 0; i < m_dim; i++){global_index[i] = i;}
	int *local_index = mkl_malloc(sizeof *local_index * (block_height), 64);
	for(int i = 0; i < block_height; i++){local_index[i] = global_index[i + rank*block_height];}

	/* Allocating memory for k-i pair tracking; one set for global tracking, and one set for local tracking. */
	int *global_k = mkl_malloc(sizeof *global_k * (log2size*size*n_dim), 64);
        int *global_i = mkl_malloc(sizeof *global_i * (log2size*size*n_dim), 64);
	int *local_k = mkl_malloc(sizeof *local_k * (n_dim), 64);
        int *local_i = mkl_malloc(sizeof *local_i * (n_dim), 64);

	/* Allocating memory for the receive buffer. */
	double *W_buffer = mkl_malloc(sizeof *W_buffer * (b*n_dim), 64);	
	double end_malloc  = MPI_Wtime();

	/*
 	 * 	BEGIN TOURNAMENT PIVOTING LOOP
 	 */ 	
	double start_loop = MPI_Wtime();
	/* Loop passing through layers of the BRT. */
	for (int d = 0; d < log2(lastpower); d++) 
	{
		/* Setting the values of the local k-i pairs to -1. */
		for(int i = 0; i < n_dim; i++){local_k[i] = -1; local_i[i] = -1;}
	
			if(d==0){
				/* Performing communication between senders and receivers in the first layer. */
				for (int k = 0; k < lastpower; k += 1 << (d + 1))
				{
					const int receiver = k;
					const int sender = k + (1 << d);
					if (rank == receiver)
					{
						MPI_Recv(W_buffer, b*n_dim, MPI_DOUBLE, sender, tag1, MPI_COMM_WORLD, &status); /* Receive "winning" rows from neighbours tournament. */
						pivot_contestants(A, block_height,n_dim, local_index, local_k, local_i); 	/* Perform pivoting between eligible entries of sub-matrix 
																 * currently located on processor. */
						double *W = collect_winners(A,b,n_dim); 					/* Collect winning rows from this processors tournament. */
						MPI_Recv(&local_index[b], b, MPI_INT, sender, tag2, MPI_COMM_WORLD, &status); 	/* Receive indices of the winning rows from neighbours tournament. */
						combine_winners(A_minor, W, W_buffer, tournament_height, n_dim); 		/* Form the new sub-matrix from the winners of the round. */
						mkl_free(W);
					}
					
					else if (rank == sender)
					{
						pivot_contestants(A, block_height,n_dim, local_index, local_k, local_i); 	/* Perform pivoting between eligible entries of sub-matrix 
																 * currently located on processor. */
						MPI_Send(A, b*n_dim, MPI_DOUBLE, receiver, tag1, MPI_COMM_WORLD); 		/* Send winning rows to neighbour. */
						MPI_Send(local_index, b, MPI_INT, receiver, tag2, MPI_COMM_WORLD); 		/* Send indices of the winning rows to neighbour. */
					}
				}
			}
			else{
				/* Performing communication between senders and receivers in a specified layer. */
				for (int k = 0; k < lastpower; k += 1 << (d + 1))
				{
					const int receiver = k;
					const int sender = k + (1 << d);
					if (rank == receiver)
					{
						MPI_Recv(W_buffer, b*n_dim, MPI_DOUBLE, sender, tag1, MPI_COMM_WORLD, &status); /* Receive "winning" rows from neighbours tournament. */
						pivot_contestants(A_minor, tournament_height,n_dim, local_index, local_k, local_i);        /* Perform pivoting between eligible entries of sub-matrix 
																 * currently located on processor. */
						double *W2 = collect_winners(A_minor,b,n_dim);                                  /* Collect winning rows from this processors tournament. */
						MPI_Recv(&local_index[b], b, MPI_INT, sender, tag2, MPI_COMM_WORLD, &status);   /* Receive indices of the winning rows from neighbours tournament. */
						combine_winners(A_minor, W2, W_buffer, tournament_height , n_dim);             	/* Form the new sub-matrix from the winners of the round. */
						mkl_free(W2);
					}

					else if (rank == sender)
					{
						pivot_contestants(A_minor, tournament_height,n_dim, local_index, local_k, local_i);        /* Perform pivoting between eligible entries of sub-matrix 
																 * currently located on processor. */
						MPI_Send(A_minor, b*n_dim, MPI_DOUBLE, receiver, tag1, MPI_COMM_WORLD);         /* Send winning rows to neighbour. */
						MPI_Send(local_index, b, MPI_INT, receiver, tag2, MPI_COMM_WORLD);              /* Send indices of the winning rows to neighbour. */
					}
				}
				
			}
		/* Barrier to ensure that all ranks have completed their round of tournament pivoting. */
		MPI_Barrier(MPI_COMM_WORLD);
		
		/* Rank 0 collects all pivoting information from a given round into the global k-i array */
		if(rank == 0)
		{
			for(int s = 1; s < size; s++)
			{
				MPI_Recv(&global_k[(s*n_dim)+(size*d*n_dim)], n_dim, MPI_INT, s, tag3, MPI_COMM_WORLD, &status);
				MPI_Recv(&global_i[(s*n_dim)+(size*d*n_dim)], n_dim, MPI_INT, s, tag4, MPI_COMM_WORLD, &status);
			}

			/* Rank 0 must also record its own pivots for a given round. */
			for(int i = 0; i < n_dim; i++)
                        {
                                global_k[(d*size*n_dim)+i] = local_k[i];
                                global_i[(d*size*n_dim)+i] = local_i[i];
                        }
		}

		/* All non-root ranks send their pivoting information to rank 0 to track pivots being used. */
		else if(rank != 0)
		{
			MPI_Send(local_k, n_dim, MPI_INT, 0, tag3, MPI_COMM_WORLD);
			MPI_Send(local_i, n_dim, MPI_INT, 0, tag4, MPI_COMM_WORLD);
		}

        }
	/* Final barrier to ensure all communications are complete at this point. */
	MPI_Barrier(MPI_COMM_WORLD);
	double end_loop  = MPI_Wtime();
	/*
 	 *	END TOURNAMENT PIVOTING LOOP
 	 */ 	
	
	double start_pivoting = MPI_Wtime();
	if(rank == 0)
	{
		/* One final round of pivoting takes place to determine the overall best pivots according to tournament pivoting. */
		if(size>1){pivot_contestants(A_minor, tournament_height, n_dim, local_index, local_k, local_i);}
		else{pivot_contestants(A, block_height, n_dim, local_index, local_k, local_i);}	
		/* The pivots that have been collected across all of the different rounds are now collected and applied to the inital TS matrix */
		
		perform_pivoting2(global_index, global_k, global_i, local_k, local_i, m_dim, n_dim, log2size, size);
		MPI_Scatter(global_index, block_height, MPI_INT, local_index, block_height, MPI_INT, rank, MPI_COMM_WORLD);	
	}

	permute_matrix(mat, local_index, A, block_height, n_dim);
	MPI_Gather(A, block_size, MPI_DOUBLE, mat, block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	double end_pivoting = MPI_Wtime();

	if(rank == 0)
	{
		/* Finally, we can peform LU factorisation on the permuted TS matrix, resulting in a successful implementation of TSLU. */	
		double start_ge = MPI_Wtime();
		gaussian_elimination(mat, L, U, m_dim, n_dim);		
		double end_ge = MPI_Wtime();
		
		double *L_new = reshape_L(L, m_dim, n_dim); /* Converts L to an [mxm] matrix to multiply against U */
                calculate_LU(LU, L_new, U, m_dim, n_dim); mkl_free(L_new);      /* Calculating the product of L and U */

		/* Printing output matrices for visualisation */
		//lu_print_matrix(mat,m_dim,n_dim);		/* Pivoted matrix mat. */
		//lu_print_matrix(L,m_dim,n_dim);		/* L from LU factorisation. */
		//lu_print_matrix(U,m_dim,n_dim);		/* U from LU factorisation. */
		//lu_print_matrix(LU,m_dim,n_dim);         	/* Recombined LU to compare to A. */
	
		double start_stability = MPI_Wtime();
		double permuted_mat_inf_norm = infinity_norm(mat, m_dim, n_dim);
		double LU_inf_norm = infinity_norm(LU, m_dim, n_dim);
		const double growth_factor = LU_inf_norm/permuted_mat_inf_norm;		
		const double backwards_stability = (permuted_mat_inf_norm - LU_inf_norm)/original_mat_inf_norm;
		double end_stability = MPI_Wtime();

		/* Collecting runtime */
                double end = MPI_Wtime();
		
		/*
		const double runtime_malloc = end_malloc-start_malloc;
		const double runtime_loop = end_loop-start_loop;
		const double runtime_pivoting = end_pivoting-start_pivoting;
		const double runtime_ge = end_ge-start_ge;
		const double runtime_stability = end_stability-start_stability;
		const double runtime = end-start;
		*/		
		const double runtime_adj = end_pivoting - start_loop;	/* The adjusted runtime times only the duration of the pivoting operations. */

		print_outputs(matrix_build_type, filename, m_dim, n_dim, size, runtime, growth_factor, backwards_stability);

	}	
		
	/* Cleaning up memory */ 
	mkl_free(mat);
	mkl_free(A);
	mkl_free(A_minor);
	mkl_free(L);
	mkl_free(U);
	mkl_free(LU);
	mkl_free(global_index);
	mkl_free(local_index);
	mkl_free(global_k); 
	mkl_free(global_i);
	mkl_free(local_k);
	mkl_free(local_i);
	mkl_free(W_buffer);

	/* Finalizing MPI, and returning 0 on successful completion. */
	MPI_Finalize();
	return 0;
}
