/**
 * @file 	file_ops.h
 * @brief 	Contains function declarations relating to opening
 *              matrices present in the matlab generated csv files.
 * @author 	William O'Sullivan
 * @version 	1.2
 * @date 	2021-09-03
 */

#ifndef FILE_OPS_H_2021_03_26_14_23
#define FILE_OPS_H_2021_03_26_14_23

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <time.h>
#include "mkl.h"

/* File read: */
void read_panel(double * const mat, const int m_dim, const int n_dim, const char * filename);

/* Create output: */
void print_outputs(const int matrix_build_type, const char * matrix_type, const int m_dim, const int n_dim, const int size, const double runtime, const double growth_factor, const double backwards_stability);

#endif /* End include guard: FILE_OPS_H_2021_03_26_14_23 */
