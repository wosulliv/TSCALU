#!/bin/bash

for file in *.zip;
do
	unzip ${file} -d matrices_4096_zip_and_tar;
done
