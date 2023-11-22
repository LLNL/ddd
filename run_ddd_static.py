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

# if any(fname.endswith('.csv') for fname in os.listdir()):
#     print("==> Current working directory is", os.getcwd())
#     print("==> Remove all .csv files from current working directory, then run again.")
#     exit()

# jobid = 123456
# zoomstep = 0.4
# minzoom = 0.0
# maxzoom = 4.0
# numtrials = 10
# zoomexps = np.linspace(minzoom,maxzoom,num=int(maxzoom/zoomstep))

# for zoomexp in zoomexps:
#     for seed in range(numtrials):
#         print("\n ==> Starting ddd trial with zoom exponent =",zoomexp, " seed=", seed, "\n")
#         subprocess.run(["python", "delaunay_density_diagnostic.py", "--jobid", str(jobid), "--seed", str(seed), "--zoomexp", str(zoomexp)])
#     #
# #

allfiles = []
for x in os.listdir():
    if x.endswith(".csv"):
        allfiles.append(str(x))

n_fnames = ["{}\n".format(i) for i in allfiles]
with open('allfiles.multi', 'w') as fp:
    fp.writelines(n_fnames)

subprocess.run(["python", "generate_ddd_figures.py", "allfiles.multi"])
