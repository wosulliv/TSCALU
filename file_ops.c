/**
 * @file 	file_ops.c
 * @brief 	Contains functions relating to opening
 * 		matrices present in the matlab generated
 * 		csv files.
 * @author 	William O'Sullivan
 * @version 	1.2
 * @date 	2021-08-30
 */

/* IMPORTING LIBRARIES/MODULES/HEADER FILES */
#include"file_ops.h"

/**
 * @brief 	Creates a panel W of size [mxn] 
 * @params[in] 	mat		Matrix memory.
 * @params[in] 	m_dim		Number of rows in the panel.
 * @params[in] 	n_dim		Number of columns in the panel.
 * @params[in] 	filename	Name of the file containing the gallery matrix of interest.
 */
void read_panel(double * const mat, const int m_dim, const int n_dim, const char * filename)
{
	FILE *fp; 			/* Defining file pointer. */
	const int true_dim = 4096; 	/* Size of square generated matrices. */
	double *full_mat = mkl_malloc(sizeof(*full_mat) * true_dim * true_dim, 64);

	/* Open "panel.txt" in write mode. */
        if((fp = fopen(filename, "r")) == NULL){perror("Failed to open file for reading."); exit(EXIT_FAILURE);}

	/* Reading file into internal object. */
        for(int i = 0; i < true_dim; i++){for(int j = 0; j < true_dim; j++){fscanf(fp, "%lf", full_mat+i*true_dim+j);}}
		
	/* Reading memory into the matrix format. */
	for(int i = 0; i < m_dim; i++){for(int j = 0; j < n_dim; j++){mat[i*n_dim+j] = full_mat[i*(m_dim)+j];}}	
	
	fclose(fp);
	mkl_free(full_mat);
}

/**
 * @brief       Formats and collects all of the relevant outputs for analysis.
 * @params[in]	matrix_build_type	Describes whether the matrix was a special variant or not.
 * @params[in]	matrix_type		Gives the name of of the special variant.
 * @params[in]	m_dim			Number of rows in the panel.
 * @params[in]	n_dim			Number of columns in the panel.
 * @params[in]	size			Number of processors used in the assessment.
 * @params[in]	runtime			Time taken to process main loop and perform pivoting and GE.
 * @params[in]	growth_factor		Measures deviation between PA and LU
 * @params[in]	backwards_stability	Analyses difference between PA and LU relative to A
 */
void print_outputs(const int matrix_build_type, const char * matrix_type, const int m_dim, const int n_dim, const int size, const double runtime, const double growth_factor, const double backwards_stability)
{
	FILE *fp;
	int tabular = 0;

	char file_name[64];
	struct tm *timenow;

	time_t now = time(NULL);
	timenow = gmtime(&now);

	strftime(file_name, sizeof(file_name), "results/outfile_%Y-%m-%d_%H:%M:%S.txt",timenow);

	if((fp = fopen(file_name, "w")) == NULL) /* open outfile_*.txt in write mode */
        {
                perror("Failed to open file for writing");
                exit(EXIT_FAILURE);
        }

	if(tabular == 0)
	{
		fprintf(fp, "Matrix Type:\n\t");
		if(matrix_build_type==0)
		{
			fprintf(fp, "Random Gaussian Matrix\n");
		}
		else
		{
			fprintf(fp, matrix_type);
			fprintf(fp, "\n");
		}

		fprintf(fp, "Panel Size [mxn]:\n\t");
		fprintf(fp, "%d  by  %d\n", m_dim, n_dim);

		fprintf(fp, "Number of Processors Used:\n\t");
		fprintf(fp, "%d procs\n", size);

		fprintf(fp, "Runtime:\n\t");
		fprintf(fp, "%f s\n", runtime);

		fprintf(fp, "Growth Factor:\n\t");
		fprintf(fp, "%f\n", growth_factor);

		fprintf(fp, "Backwards Stability:\n\t");
                fprintf(fp, "%f\n", backwards_stability);
	}

	if(tabular == 1)
	{
		if(matrix_build_type==0)
                {
                        fprintf(fp, "Random Gaussian Matrix");
                }
                else
                {
                        fprintf(fp, matrix_type);
                }
		fprintf(fp, "\n%d\n", m_dim);
		fprintf(fp, "%d\n", n_dim);
		fprintf(fp, "%d\n", size);
		fprintf(fp, "%f\n", runtime);
		fprintf(fp, "%f\n", growth_factor);
		//fprintf(fp, "%f\n", backwards_stability);
	}
	fclose(fp);
}

