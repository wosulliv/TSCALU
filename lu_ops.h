/**
 * @file 	lu_ops.h
 * @brief 	Contains function declarations relating to TS 
 * 		matrix population and manipulation, as well
 * 		as the suite of functions required for 
 * 		tournament pivoting. Additionally, this
 * 		header file also prototypes a function for
 * 		explicit LUP factorisation.
 * @author 	William O'Sullivan
 * @version 	1.9
 * @date 	2021-07-14
 */

#ifndef LU_OPS_H_2021_03_26_14_23
#define LU_OPS_H_2021_03_26_14_23

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "mkl.h"

/* Matrix population and printing: */
int lu_build_matrix(double * const matrix_memory, const int m_dim, const int n_dim);
int lu_print_matrix(double * const matrix_memory, const int m_dim, const int n_dim);

/* Explicit LUP factorisation: */
void LUP(double * const A, double * const L, double * const U, double * const P, const int m_dim, const int n_dim);

/* Functions for tournament pivoting: */
void pivot_contestants(double * const A, const int m_dim, const int n_dim, int * const local_index, int * const local_k, int * const local_i);
double * collect_winners(double * const A, const int m_dim, const int n_dim);
void combine_winners(double * const A, double * const W, double * const W_buff, const int m_dim, const int n_dim);

void perform_pivoting(double * const mat, int * const global_index, int * const global_k, int * const global_i, int * const local_k, int * const local_i, const int m_dim, const int n_dim, const int log2size, const int size);
void swap_indices(int * const global_index, int * const global_index_tmp, int k, int i);

/* Function for performing LU factorisation on the pivoted TS matrix A: */
void gaussian_elimination(double * const A, double * const L, double * const U, const int m_dim, const int n_dim);

/* Function for reshaping the matrix U: */
double * reshape_U(double * const U, const int n_dim);

/* Function for reshaping the matrix L: */
double * reshape_L(double * const L, const int m_dim, const int n_dim);

/* Functions for matrix multipliction: */
void calculate_PA(double * const PA, double * const P, double * const A, const int m_dim, const int n_dim);
void calculate_LU(double * const LU, double * const L, double * const U, const int m_dim, const int n_dim);

/* Functions for growth factor determination: */
double infinity_norm(double * const matrix, const int m_dim, const int n_dim);


#endif /* End include guard: LU_OPS_H_2021_03_26_14_23 */
