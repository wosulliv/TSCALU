# -*- coding: utf-8 -*-
###############################################################################
#University:        Trinity College Dublin
#Course:            High Performance Computing (MSc)
#Assignment:        Dissertation
#Supervisor:        Kirk Soodhalter
#Author:            William O'Sullivan
#Student No.:       16321101
#Contact email:     wosulliv@tcd.ie
#Created:           2021-09-10
#Due date:          2021-09-17
#IDE Used:          Spyder 4.1.5
#Version:           Python 3.7
###############################################################################

import csv
import matplotlib.pyplot as plt

#txt_title = file_name + ".txt"
with open("strong.csv") as file:
    parsing = csv.reader(file)
    data = list(parsing)
data[0][0] = 'name'            
    
gaussian_core = []
gaussian_time = []
circulant_core = []
circulant_time = []
diag_dom_core = []
diag_dom_time = []
jord_block_core = []
jord_block_time = []
kms_core = []
kms_time = []
neumann_core = []
neumann_time = []

for entry in data:
    if entry[0] == "Gaussian":
        gaussian_core.append(float(entry[3]))
        gaussian_time.append(float(entry[4]))
    if entry[0] == "Circulant":
        circulant_core.append(float(entry[3]))
        circulant_time.append(float(entry[4]))
    if entry[0] == "Diagonally Dominant":
        diag_dom_core.append(float(entry[3]))
        diag_dom_time.append(float(entry[4]))
    if entry[0] == "Jordan Block":
        jord_block_core.append(float(entry[3]))
        jord_block_time.append(float(entry[4]))
    if entry[0] == "KMS":
        kms_core.append(float(entry[3]))
        kms_time.append(float(entry[4]))
    if entry[0] == "Neumann":
        neumann_core.append(float(entry[3]))
        neumann_time.append(float(entry[4]))
 
gaussian_speedup = []
circulant_speedup = []
diag_dom_speedup = []
jord_block_speedup = []
kms_speedup = []
neumann_speedup = []

for i in range(4):
    gaussian_speedup.append((gaussian_time[0]/gaussian_time[i]))
    circulant_speedup.append((circulant_time[0]/circulant_time[i]))
    diag_dom_speedup.append((diag_dom_time[0]/diag_dom_time[i]))
    jord_block_speedup.append((jord_block_time[0]/jord_block_time[i]))
    kms_speedup.append((kms_time[0]/kms_time[i]))
    neumann_speedup.append((neumann_time[0]/neumann_time[i]))
    
ideal_core = []
ideal_time = []
errorbar_y = []
for i in range(4):
    ideal_core.append(gaussian_core[i])
    ideal_time.append(gaussian_core[i])
    errorbar_y.append(0.1)


    
fig_1 = plt.figure(figsize=(12,6.75))
ax_1 = fig_1.add_subplot(1,1,1)
title_text_1 = "Strong Scaling in Different Matrices [4096x512]"
ax_1.set_title(title_text_1, fontweight="bold", fontsize=15)
ax_1.set_xlabel('Cores', fontweight="bold", fontsize=15)
ax_1.set_ylabel('Speedup', fontweight="bold", fontsize=15)
#ax_1.set_xlim([0.9,8.1])
#ax_1.set_ylim([0.9,8.1])
ax_1.plot(gaussian_core, gaussian_speedup, 'r:')
ax_1.plot(gaussian_core, gaussian_speedup, 'ro', label='Random Gaussian Matrix Speedup')
ax_1.plot(circulant_core, circulant_speedup, 'g--')
ax_1.plot(circulant_core, circulant_speedup, 'gs', label='Circulant Matrix Speedup')
ax_1.plot(diag_dom_core, diag_dom_speedup, 'b-.')
ax_1.plot(diag_dom_core, diag_dom_speedup, 'b^', label='Diagonally Dominant Matrix Speedup')
ax_1.plot(jord_block_core, jord_block_speedup, 'c:')
ax_1.plot(jord_block_core, jord_block_speedup, 'cD', label='Jordan Block Matrix Speedup')
ax_1.plot(kms_core, kms_speedup, 'm--')
ax_1.plot(kms_core, kms_speedup, 'mX', label='KMS Matrix Speedup')
# ax_1.plot(kms_core, kms_speedup, 'm--')
# ax_1.errorbar(kms_core, kms_speedup, yerr=errorbar_y, capsize=5, linestyle='-', label='KMS Matrix Speedup')
ax_1.plot(neumann_core, neumann_speedup, 'y-.')
ax_1.plot(neumann_core, neumann_speedup, 'yP', label='Neumann Matrix Speedup')
#ax_1.plot(ideal_core, ideal_time, 'k*-', label='Ideal Speedup')
ax_1.legend(fontsize=15, loc='best',fancybox= True,framealpha=0.80)
ax_1.grid()

#%%
with open("weak.csv") as file:
    parsing = csv.reader(file)
    data = list(parsing)
data[0][0] = 'name'            
    
gaussian_core = []
gaussian_time = []
circulant_core = []
circulant_time = []
diag_dom_core = []
diag_dom_time = []
jord_block_core = []
jord_block_time = []
kms_core = []
kms_time = []
neumann_core = []
neumann_time = []

for entry in data:
    if entry[0] == "Gaussian":
        gaussian_core.append(float(entry[3]))
        gaussian_time.append(float(entry[4]))
    if entry[0] == "Circulant":
        circulant_core.append(float(entry[3]))
        circulant_time.append(float(entry[4]))
    if entry[0] == "Diagonally Dominant":
        diag_dom_core.append(float(entry[3]))
        diag_dom_time.append(float(entry[4]))
    if entry[0] == "Jordan Block":
        jord_block_core.append(float(entry[3]))
        jord_block_time.append(float(entry[4]))
    if entry[0] == "KMS":
        kms_core.append(float(entry[3]))
        kms_time.append(float(entry[4]))
    if entry[0] == "Neumann":
        neumann_core.append(float(entry[3]))
        neumann_time.append(float(entry[4]))
 
gaussian_speedup = []
circulant_speedup = []
diag_dom_speedup = []
jord_block_speedup = []
kms_speedup = []
neumann_speedup = []

for i in range(4):
    gaussian_speedup.append((gaussian_time[0]/gaussian_time[i])/gaussian_core[i])
    circulant_speedup.append((circulant_time[0]/circulant_time[i])/gaussian_core[i])
    diag_dom_speedup.append((diag_dom_time[0]/diag_dom_time[i])/gaussian_core[i])
    jord_block_speedup.append((jord_block_time[0]/jord_block_time[i])/gaussian_core[i])
    kms_speedup.append((kms_time[0]/kms_time[i])/gaussian_core[i])
    neumann_speedup.append((neumann_time[0]/neumann_time[i])/gaussian_core[i])
    
ideal_core = []
ideal_time = []
for i in range(4):
    ideal_core.append(gaussian_core[i])
    ideal_time.append(gaussian_core[i])

    
fig_2 = plt.figure(figsize=(12,6.75))
ax_1 = fig_2.add_subplot(1,1,1)
title_text_1 = "Weak Scaling in Different Matrices [512x512] at per core"
ax_1.set_title(title_text_1, fontweight="bold", fontsize=15)
ax_1.set_xlabel('Cores', fontweight="bold", fontsize=15)
ax_1.set_ylabel('Speedup', fontweight="bold", fontsize=15)
# ax_1.set_xlim([0.9,8.1])
# ax_1.set_ylim([0.9,8.1])
ax_1.plot(gaussian_core, gaussian_speedup, 'r:')
ax_1.plot(gaussian_core, gaussian_speedup, 'ro', label='Random Gaussian Matrix Speedup')
ax_1.plot(circulant_core, circulant_speedup, 'g--')
ax_1.plot(circulant_core, circulant_speedup, 'gs', label='Circulant Matrix Speedup')
ax_1.plot(diag_dom_core, diag_dom_speedup, 'b-.')
ax_1.plot(diag_dom_core, diag_dom_speedup, 'b^', label='Diagonally Dominant Matrix Speedup')
ax_1.plot(jord_block_core, jord_block_speedup, 'c:')
ax_1.plot(jord_block_core, jord_block_speedup, 'cD', label='Jordan Block Matrix Speedup')
ax_1.plot(kms_core, kms_speedup, 'm--')
ax_1.plot(kms_core, kms_speedup, 'mX', label='KMS Matrix Speedup')
ax_1.plot(neumann_core, neumann_speedup, 'y-.')
ax_1.plot(neumann_core, neumann_speedup, 'yP', label='Neumann Matrix Speedup')
#ax_1.plot(ideal_core, ideal_time, 'k*-', label='Ideal Speedup')
ax_1.legend(fontsize=15, loc='best',fancybox= True,framealpha=0.80)
ax_1.grid()

#%%
with open("RGM_strong.csv") as file:
    parsing = csv.reader(file)
    data = list(parsing)
data[0][0] = 'name'            
    
gaussian_core = []
gaussian_time = []


for entry in data:
    if entry[0] == "Gaussian":
        gaussian_core.append(float(entry[3]))
        gaussian_time.append(float(entry[4]))
 
gaussian_speedup = []


for i in range(4):
    gaussian_speedup.append((gaussian_time[0]/gaussian_time[i]))
    
ideal_core = []
ideal_time = []
for i in range(4):
    ideal_core.append(gaussian_core[i])
    ideal_time.append(gaussian_core[i])

    
fig_3 = plt.figure(figsize=(12,6.75))
ax_1 = fig_3.add_subplot(1,1,1)
title_text_1 = "Strong Scaling in Random Gaussian Matrices [16384x2048]"
ax_1.set_title(title_text_1, fontweight="bold", fontsize=15)
ax_1.set_xlabel('Cores', fontweight="bold", fontsize=15)
ax_1.set_ylabel('Speedup', fontweight="bold", fontsize=15)
#ax_1.set_xlim([0.9,8.1])
#ax_1.set_ylim([0.9,8.1])
ax_1.plot(gaussian_core, gaussian_speedup, 'r:')
ax_1.plot(gaussian_core, gaussian_speedup, 'ro', label='Random Gaussian Matrix Speedup')
# ax_1.plot(circulant_core, circulant_speedup, 'g--')
# ax_1.plot(circulant_core, circulant_speedup, 'gs', label='Circulant Matrix Speedup')
# ax_1.plot(diag_dom_core, diag_dom_speedup, 'b-.')
# ax_1.plot(diag_dom_core, diag_dom_speedup, 'b^', label='Diagonally Dominant Matrix Speedup')
# ax_1.plot(jord_block_core, jord_block_speedup, 'c:')
# ax_1.plot(jord_block_core, jord_block_speedup, 'cD', label='Jordan Block Matrix Speedup')
# ax_1.plot(kms_core, kms_speedup, 'm--')
# ax_1.plot(kms_core, kms_speedup, 'mX', label='KMS Matrix Speedup')
# ax_1.plot(neumann_core, neumann_speedup, 'y-.')
# ax_1.plot(neumann_core, neumann_speedup, 'yP', label='Neumann Matrix Speedup')
#ax_1.plot(ideal_core, ideal_time, 'k*-', label='Ideal Speedup')
ax_1.legend(fontsize=15, loc='best',fancybox= True,framealpha=0.80)
ax_1.grid()

#%%
with open("RGM_weak.csv") as file:
    parsing = csv.reader(file)
    data = list(parsing)
data[0][0] = 'name'            
    
gaussian_core = []
gaussian_time = []


for entry in data:
    if entry[0] == "Gaussian":
        gaussian_core.append(float(entry[3]))
        gaussian_time.append(float(entry[4]))
 
gaussian_speedup = []


for i in range(4):
    gaussian_speedup.append((gaussian_time[0]/gaussian_time[i])/gaussian_core[i])
    
ideal_core = []
ideal_time = []
for i in range(4):
    ideal_core.append(gaussian_core[i])
    ideal_time.append(gaussian_core[i])

    
fig_4 = plt.figure(figsize=(12,6.75))
ax_1 = fig_4.add_subplot(1,1,1)
title_text_1 = "Weak Scaling in Random Gaussian Matrices [2048x2048] per core"
ax_1.set_title(title_text_1, fontweight="bold", fontsize=15)
ax_1.set_xlabel('Cores', fontweight="bold", fontsize=15)
ax_1.set_ylabel('Speedup', fontweight="bold", fontsize=15)
#ax_1.set_xlim([0.9,8.1])
#ax_1.set_ylim([0.9,8.1])
ax_1.plot(gaussian_core, gaussian_speedup, 'r:')
ax_1.plot(gaussian_core, gaussian_speedup, 'ro', label='Random Gaussian Matrix Speedup')
# ax_1.plot(circulant_core, circulant_speedup, 'g--')
# ax_1.plot(circulant_core, circulant_speedup, 'gs', label='Circulant Matrix Speedup')
# ax_1.plot(diag_dom_core, diag_dom_speedup, 'b-.')
# ax_1.plot(diag_dom_core, diag_dom_speedup, 'b^', label='Diagonally Dominant Matrix Speedup')
# ax_1.plot(jord_block_core, jord_block_speedup, 'c:')
# ax_1.plot(jord_block_core, jord_block_speedup, 'cD', label='Jordan Block Matrix Speedup')
# ax_1.plot(kms_core, kms_speedup, 'm--')
# ax_1.plot(kms_core, kms_speedup, 'mX', label='KMS Matrix Speedup')
# ax_1.plot(neumann_core, neumann_speedup, 'y-.')
# ax_1.plot(neumann_core, neumann_speedup, 'yP', label='Neumann Matrix Speedup')
#ax_1.plot(ideal_core, ideal_time, 'k*-', label='Ideal Speedup')
ax_1.legend(fontsize=15, loc='best',fancybox= True,framealpha=0.80)
ax_1.grid()