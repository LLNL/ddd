#############################################################################################
# Driver script for Delaunay Density Diagnositic with static data set#
#   as described in the paper
#   Algorithm XXXX: The Delaunay Density Diagnostic
#       under review at ACM Transactions on Mathematical Software
#       (original title: Data-driven geometric scale detection via Delaunay interpolation)
#   by Andrew Gillette and Eugene Kur
#   Version 2.0, November 2023
#
# This script does the following:
#   1) Run multiple trials of delaunay_density_diagnostic.py for a static data set.
#           Output consists of one file per trial of the form zz*seed*.csv
#   2) Save list of all output files into a txt file called allfiles.multi.
#   3) Call generate_ddd_figures.py on allfiles.multi, which does the following: 
#           Generate figure showing MSD and grad-MSD rate as a function 
#               of average sample spacing.  
#           Output figure is displayed and then saved as ddd-figure.png
#############################################################################################

import subprocess
import os
import numpy as np

if any(fname.endswith('.csv') for fname in os.listdir()):
    print("==> Current working directory is", os.getcwd())
    print("==> Remove all .csv files from current working directory, then run again.")
    exit()

if any(fname.endswith('allfiles.multi') for fname in os.listdir()):
    print("==> Current working directory is", os.getcwd())
    print("==> Remove allfiles.multi from current working directory, then run again.")
    exit()

jobid = 123456
numtrials = 10
# staticdatapath = 'staticdata/examples/data_paraboloid_exp2_10000.csv'
# staticdatapath = 'staticdata/examples/data_griewank_exp2_10000.csv'
staticdatapath = 'staticdata/examples/data_griewank_dim5_10000.csv'

for seed in range(numtrials):
    print("\n ==> Starting ddd for static data =",staticdatapath, " seed=", seed+1, "\n")
    subprocess.run(["python", "delaunay_density_diagnostic.py", "--jobid", str(jobid),
                    "--staticdatapath", staticdatapath, "--seed", str(seed+1), "--numtestperdim", str(2)])
#

allfiles = []
for x in os.listdir():
    if x.endswith(".csv"):
        allfiles.append(str(x))

n_fnames = ["{}\n".format(i) for i in allfiles]
with open('allfiles.multi', 'w') as fp:
    fp.writelines(n_fnames)

subprocess.run(["python", "generate_ddd_figures.py", "allfiles.multi"])
