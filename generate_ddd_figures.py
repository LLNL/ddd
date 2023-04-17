
from matplotlib import pyplot as plt
from matplotlib.ticker import FormatStrFormatter
import math
import numpy as np
import pandas as pd
from matplotlib import colors
from matplotlib.pyplot import cm
import subprocess


import sys
try:
    # case = int(sys.argv[1])
    infile = str(sys.argv[1])
except IndexError as err:
    print("Not enough arguments: {0}".format(err))
    sys.exit(1)
except ValueError as err:
    print("Illegal value in argument: {0}".format(err))
    sys.exit(1)

min_density = 50

from functools import reduce

plt.rcParams["font.family"] = "Times New Roman"
plt.rcParams['font.size'] = 20
plt.rc('text', usetex=True)

allfiles = open(infile, 'r')
allfiles = allfiles.readlines()  
df = pd.read_csv(allfiles[0].strip(), header=0, index_col=0)   
for i in range(1,len(allfiles)):
    nextdf = pd.read_csv(allfiles[i].strip(), header=0, index_col=0) 
    df = pd.concat([df, nextdf], ignore_index=True) 

# samp0 = df['samples'].unique()[0]
# samp1 = df['samples'].unique()[1]
# samp2 = df['samples'].unique()[2]
# tempdf = df[df['samples'] == samp2]['iterate rate']
# print("val range = ", tempdf.max() - tempdf.min())
# tempdf = df[df['samples'] == samp2]['grad rate']
# print("grd range = ", tempdf.max() - tempdf.min())



###### filter out small densitites
# print("*** Filtering out densities < ", min_density)
df = df[df['density'] > min_density]
df = df.reset_index()
# print("Available sample sizes: ", df['samples'].unique())
# print("Available densities:    ", df['density'].unique())

# slurmid  = str(allfiles[0][3:9])
funcname = str(df['function name'][0])  
dimen    = str(df['dim of intput'][0])
logbase  = str(df['log base'][0])
# lb       = str(df['left bound'][0])
# rb       = str(df['right bound'][0])
tbscale  = str(np.round(df['test grid scale'][0],2))

fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(6,12))
fig.suptitle('Rates for '+funcname+', dim='+dimen+', b='+logbase, fontsize=20)

df['charL'] = df['right bound'] - df['left bound']
df['avg samp spacing'] = df['charL']/(np.power(df['samples'],1.0/float(dimen))) 
# print("Available spacings:    ", df['avg samp spacing'].unique())

###### compute zoom breaks
md = df['density'].unique().max()
zoom_breaks = df[df['density'] == md]['avg samp spacing'].unique()


for i in range(2):
    if i == 0:
        rate_to_plot = 'iterate rate'
    else: 
        rate_to_plot = 'grad rate'        

    dfsub = df[[rate_to_plot, 'avg samp spacing']]

    ir_mean   =dfsub.groupby('avg samp spacing').mean().rename(columns={"iterate rate norm inside": "ir_mean"})
    ir_25per  =dfsub.groupby('avg samp spacing').quantile(q=.25).rename(columns={"iterate rate norm inside": "ir_25per"})
    ir_75per  =dfsub.groupby('avg samp spacing').quantile(q=.75).rename(columns={"iterate rate norm inside": "ir_75per"})
    ir_max    =dfsub.groupby('avg samp spacing').max().rename(columns={"iterate rate norm inside": "ir_max"})
    ir_min    =dfsub.groupby('avg samp spacing').min().rename(columns={"iterate rate norm inside": "ir_min"})
    ir_10per  =dfsub.groupby('avg samp spacing').quantile(q=.10).rename(columns={"iterate rate norm inside": "ir_10per"})
    ir_90per  =dfsub.groupby('avg samp spacing').quantile(q=.90).rename(columns={"iterate rate norm inside": "ir_90per"})
    # ir_5per   =dfsub.groupby('avg samp spacing').quantile(q=.05).rename(columns={"iterate rate norm inside": "ir_5per"})
    # ir_95per  =dfsub.groupby('avg samp spacing').quantile(q=.95).rename(columns={"iterate rate norm inside": "ir_95per"})
    
    groupP = dfsub.groupby('avg samp spacing')
    #groupby attributes
    ir_mean        =groupP.mean().rename(columns={rate_to_plot: "ir_mean"})
    ir_25per       =groupP.quantile(q=.25).rename(columns={rate_to_plot: "ir_25per"})
    ir_75per       =groupP.quantile(q=.75).rename(columns={rate_to_plot: "ir_75per"})
    ir_max         =groupP.max().rename(columns={rate_to_plot: "ir_max"})
    ir_min         =groupP.min().rename(columns={rate_to_plot: "ir_min"})
    ir_10per       =groupP.quantile(q=.10).rename(columns={rate_to_plot: "ir_10per"})
    ir_90per       =groupP.quantile(q=.90).rename(columns={rate_to_plot: "ir_90per"})
    # ir_5per        =groupP.quantile(q=.05).rename(columns={rate_to_plot: "ir_5per"})
    # ir_95per       =groupP.quantile(q=.95).rename(columns={rate_to_plot: "ir_95per"})
    by_sample_count=groupP.count().rename(columns={rate_to_plot:"by_sample_count"})

    # Make a list of the dataframes
    data_frames = [ir_mean, ir_25per, ir_75per, ir_max, ir_min, 
                ir_10per, ir_90per, by_sample_count]

    # Merge them all at once
    merged_df = pd.concat(data_frames, join='outer', axis=1)

    x     = merged_df.index

    y_90  = np.ma.masked_values(merged_df.ir_90per, 0)
    y_75  = np.ma.masked_values(merged_df.ir_75per, 0)
    y     = np.ma.masked_values(merged_df.ir_mean , 0)
    y_25  = np.ma.masked_values(merged_df.ir_25per, 0)
    y_10  = np.ma.masked_values(merged_df.ir_10per, 0)

    
    if (rate_to_plot == 'grad rate'):
        target     =  1*np.ones_like(x)
        noise_line = -1*np.ones_like(x)
    else:
        target     = 2*np.ones_like(x)
        noise_line = 0*np.ones_like(x)

    l_10    = ax[i].fill_between(x,y_10,y_25,color='lightblue')
    l_25    = ax[i].fill_between(x,y_25,y,color='royalblue')
    l_mean, = ax[i].plot(x,y,'k',marker='o')
    l_75    = ax[i].fill_between(x,y,y_75,color='royalblue')
    l_90    = ax[i].fill_between(x,y_75,y_90,color='lightblue')
    l_target, = ax[i].plot(x, target, '--', linewidth=4, color='tab:green')
    l_noise,  = ax[i].plot(x, noise_line, '--', linewidth=4, color='tab:red')
    # ax[i].legend([l_target,l_mean,l_75,l_90],['target rate','mean','inter-quartile range','inter-decile range'],loc=3)
    # ax[i].legend([l_target,l_mean,l_75,l_noise],['target rate','mean','inter-quartile range','noise-only rate'],prop={'size': 12},loc=3)
    # ax[i].legend([l_target,l_noise,l_mean,l_75,l_90],['recoverable features','noisy features','mean rate','inter-quartile range','inter-decile range'],loc=7)
    ax[i].set_xscale('log')
    
    # draw zoom breaks
    for zb in zoom_breaks:
        ax[i].axvline(x=zb)  

    if i == 0:
        ax[i].set_ylabel(r'\texttt{MSD} rate', fontsize=24) 
        ax[i].set_ylim(-1.0, 3.0)
        ax[i].set_yticks([-1,0,1,2,3])
    else:
        ax[i].set_ylabel(r'\texttt{grad-MSD} rate', fontsize=24)
        ax[i].set_ylim(-2.0, 2.0)
        ax[i].set_yticks([-2,-1,0,1,2])
        ax[i].set_xlabel('average sample spacing', fontsize=24) # L/N^(1/d)
        
    plt.xticks(fontsize=20)
    plt.yticks(fontsize=20)



plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.05)

figname = 'ddd-figure.png'
plt.savefig(figname,format='png', bbox_inches='tight')
# plt.show()
print("==> Saved figure as ddd-figure.png")
subprocess.call(["open", figname])