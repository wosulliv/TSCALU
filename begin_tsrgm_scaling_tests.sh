#!/bin/bash

make;

echo "Testing RGM 16384x128 strong scaling..."

mpirun -n 1 ./run_tscalu -r 16384 -c 128 -t 0
mpirun -n 2 ./run_tscalu -r 16384 -c 128 -t 0
mpirun -n 4 ./run_tscalu -r 16384 -c 128 -t 0
mpirun -n 8 ./run_tscalu -r 16384 -c 128 -t 0

echo "Testing RGM (N*2048)x128 weak scaling..."

mpirun -n 1 ./run_tscalu -r 2048 -c 128 -t 0
mpirun -n 2 ./run_tscalu -r 4192 -c 128 -t 0
mpirun -n 4 ./run_tscalu -r 8192 -c 128 -t 0
mpirun -n 8 ./run_tscalu -r 16384 -c 128 -t 0

echo "Finished tests!"
