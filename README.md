This project has been designed to run on chuck.tchpc.tcd.ie and
requires dependencies as specified in the Makefile which you
may have to provide if testing the code elsewhere.

This project has the following dependencies:
1) cports
2) OpenMPI
3) gcc/9.2.0-gnu
4) intel/18.0.4
5) gsl

To begin running this project, start by running the command:
	./submit_job.sh
This will unzip all of the relevant matrices, run a test
assessment, and the print the results of the assessment to
the terminal.

After this, it is possible to run test assessments by running
any of the scripts that start with "begin", e.g. 
	./begin_sm_tests.sh

It is important to note that if you want to test any of the
special matrices, you are required to modify the 
"tscalu_main.c" script, from line 57 in order to select the
desired special matrix.

You are always able to test the Random Gaussian Matrix by
using the "-t" flag, followed by a "0".


The auxiliary scripts directory contains the Matlab scripts
used to generate and plot the special matrices, and the Python
scripts for converting the Matlab matrix .csv files into 
C-readable arrays, as well as the Python script for performance
analysis.
