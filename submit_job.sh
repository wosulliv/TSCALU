#!/bin/bash
# LOAD MODULES
module load gcc openmpi intel cports

# UNPACK MATRICES
myarray=(`find ./matrices_4096_zip_and_tar -maxdepth 1 -name "*.zip"`)
if [ ${#myarray[@]} -gt 0 ]; then 
    echo "Matrices already unzipped."
else 
    echo "Unzipping matrices."
    make matrices
fi

# RUN PROGRAM
make all
make assessment
cat results/outfile*
