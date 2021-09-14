#!/bin/bash

make;

#echo "Testing 4096x256 strong scaling..."
#
#mpirun -n 1 ./run_tscalu -r 4096 -c 256 -t 1
#mpirun -n 2 ./run_tscalu -r 4096 -c 256 -t 1
#mpirun -n 4 ./run_tscalu -r 4096 -c 256 -t 1
#mpirun -n 8 ./run_tscalu -r 4096 -c 256 -t 1
#
#echo "Testing 4096x256 weak scaling..."
#
#mpirun -n 1 ./run_tscalu -r 512 -c 256 -t 1
#mpirun -n 2 ./run_tscalu -r 1024 -c 256 -t 1
#mpirun -n 4 ./run_tscalu -r 2048 -c 256 -t 1
#mpirun -n 8 ./run_tscalu -r 4096 -c 256 -t 1
#
#echo "Finished tests!"

echo "Testing rgm 4096x256 strong scaling..."

mpirun -n 1 ./run_tscalu -r 4096 -c 256 -t 0
mpirun -n 2 ./run_tscalu -r 4096 -c 256 -t 0
mpirun -n 4 ./run_tscalu -r 4096 -c 256 -t 0
mpirun -n 8 ./run_tscalu -r 4096 -c 256 -t 0

echo "Testing rgm 4096x256 weak scaling..."

mpirun -n 1 ./run_tscalu -r 512 -c 256 -t 0
mpirun -n 2 ./run_tscalu -r 1024 -c 256 -t 0
mpirun -n 4 ./run_tscalu -r 2048 -c 256 -t 0
mpirun -n 8 ./run_tscalu -r 4096 -c 256 -t 0

echo "Finished tests!"
