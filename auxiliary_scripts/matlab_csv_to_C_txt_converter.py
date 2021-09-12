# -*- coding: utf-8 -*-
###############################################################################
#University:        Trinity College Dublin
#Course:            High Performance Computing (MSc)
#Assignment:        Dissertation
#Supervisor:        Kirk Soodhalter
#Author:            William O'Sullivan
#Student No.:       16321101
#Contact email:     wosulliv@tcd.ie
#Created:           2021-09-03
#Due date:          2021-09-17
#IDE Used:          Spyder 4.1.5
#Version:           Python 3.7
###############################################################################
"""
To convert matlab.csv matrices into arrays parsable by C
"""
import os
import csv

file_list = os.listdir()
file_name_list = []
for file_name in file_list:
    if '.csv' in file_name:
        file_name_list.append(file_name[:-4])

###
def convert_csv_to_txt(file_name):
    csv_title = file_name + ".csv"
    with open(csv_title) as file:
        parsing = csv.reader(file)
        data = list(parsing)
    
    flat_list = [item for sublist in data for item in sublist]
    float_list = []
    for i in range(0,len(flat_list)):
        float_list.append(float(flat_list[i]))

    txt_title = file_name + ".txt"
    file_pointer = open(txt_title, "w")
    for val in float_list:
        file_pointer.write(str(val))
        file_pointer.write("\n")
    file_pointer.close()
###

for file_name in file_name_list:
    convert_csv_to_txt(file_name)