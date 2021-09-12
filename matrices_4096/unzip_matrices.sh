#!/bin/bash

for file in matrices_4096_zip_and_tar/*.zip;
do
	unzip ${file} -d matrices_4096;
done
