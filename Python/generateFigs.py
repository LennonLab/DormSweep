from __future__ import division
import pandas as pd
import  matplotlib.pyplot as plt
import os, argparse

mydir = os.path.expanduser('~/github/DormSweep/')

def sweepFig():
    #Ss = [0.001, 0.01, 0.1]
    Ss = [0.1]
    for s in Ss:
        fig, ax = plt.subplots()
        g = []
        p = []
        df = pd.read_csv(mydir + 'data/sweep_' + str(s) + '.txt', header = None)
        #diff_g = []
        #diff_
        for index, row in df.iterrows():
            g.extend(row[~row.isnull()].index.values)
            p.extend(row[~row.isnull()].values)
        #for x in range(0, max(g)):
        ax.scatter(g, p, lw=2, color='blue', alpha = 0.1)
        ax.set_xlim([0,max(g)])
        ax.set_ylim([0,1])
        ax.set_xlabel('Generation', fontsize=20)
        ax.set_ylabel('Frequency of favored allele', fontsize=20)
        fig.savefig(mydir + 'figs/test.png', \
            bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
        plt.close()

sweepFig()
