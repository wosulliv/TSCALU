CC=mpicc

CFLAGS= -Wextra -Wall -g

LDFLAGS= -I{$(MKLROOT)/include} -L{$(MKLROOT)/lib/intel64} -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lpthread -lgsl -lgslcblas -lm -llapack

EXECS=  run_tscalu
MATRICES= unzip_matrix_directory unzip_matrices

all: $(EXECS)

matrices: $(MATRICES)

run_tscalu: tscalu_main.c lu_ops.o file_ops.o
	${CC} -o $@ $^ $(CFLAGS) $(LDFLAGS)

lu_ops.o: lu_ops.c
	${CC} -c $< $(CFLAGS) $(LDFLAGS)

file_ops.o: file_ops.c
	${CC} -c $< $(CFLAGS) $(LDFLAGS)

unzip_matrices: matrices_4096/unzip_matrices.sh matrices_4096_zip_and_tar/unzip_matrix_directory.sh
	matrices_4096/unzip_matrices.sh

unzip_matrix_directory: matrices_4096_zip_and_tar/unzip_matrix_directory.sh
	matrices_4096_zip_and_tar/unzip_matrix_directory.sh

pre_compiled: lu_ops.h file_ops.h
	${CC} -pch -c $^ $(CFLAGS) $(LDFLAGS)

.PHONY: assessment clean

assessment:
	mpirun -n 4 ./run_tscalu -r 4096 -c 128 -t 0

clean:
	$(RM) *.o *.gch $(EXECS) matrices_4096/*.txt matrices_4096_zip_and_tar/*.zip results/outfile*.txt

