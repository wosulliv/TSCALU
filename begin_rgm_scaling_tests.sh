#!/bin/bash

make;

echo "Testing RGM 16384x1024 strong scaling..."

mpirun -n 1 ./run_tscalu -r 16384 -c 1024 -t 0
mpirun -n 2 ./run_tscalu -r 16384 -c 1024 -t 0
mpirun -n 4 ./run_tscalu -r 16384 -c 1024 -t 0
mpirun -n 8 ./run_tscalu -r 16384 -c 1024 -t 0

echo "Testing RGM (N*2048)x1024 weak scaling..."

mpirun -n 1 ./run_tscalu -r 2048 -c 1024 -t 0
mpirun -n 2 ./run_tscalu -r 4192 -c 1024 -t 0
mpirun -n 4 ./run_tscalu -r 8192 -c 1024 -t 0
mpirun -n 8 ./run_tscalu -r 16384 -c 1024 -t 0

echo "Finished tests!"
