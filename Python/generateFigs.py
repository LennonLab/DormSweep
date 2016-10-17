from __future__ import division
import pandas as pd
import  matplotlib.pyplot as plt
import os, argparse, math, itertools
import numpy as np

mydir = os.path.expanduser('~/github/DormSweep/')

def sweepFig():
    #Ss = [0.001, 0.01, 0.1]
    Ss = [0.1]
    maxg = 0
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
        p_theor = []
        for x in range(1, max(g)):
            p_theor_x = (1/1000) / ((1/1000) + ((1 - (1/1000)) * math.exp(-s*x) ) )
            p_theor.append(p_theor_x)
        ax.plot(range(1, max(g)), p_theor, lw = 2, color = 'black')
        ax.scatter(g, p, lw=2, color='blue', alpha = 0.1)

        if max(g) > maxg:
            maxg = max(g)

    ax.set_xlim([0, maxg])
    ax.set_ylim([0, 1])
    ax.set_xlabel('Generation', fontsize=20)
    ax.set_ylabel('Frequency of   allele', fontsize=20)
    fig.savefig(mydir + 'figs/test.png', \
        bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def sweepDormFig(N = 1000, M = 1000, s = 0.1):
    Cs = [1, 10, 100]
    #Cs = [100]
    maxg = 0
    fig, ax = plt.subplots()
    colors = ['#87CEEB', '#FFA500', '#FF6347']
    for j, c in enumerate(Cs):
        g = []
        p = []
        df = pd.read_csv(mydir + 'data/sweep_N_' + str(N) + '_M_' + str(M) + '_c_' + \
            str(c) + '_s_' + str(s) + '.txt', header = None)

        for index, row in df.iterrows():
            g.extend(row[~row.isnull()].index.values)
            p.append(row[~row.isnull()].values)

        #g_plot = [np.mean(x) for x in zip(*g)]
        #print itertools.izip_longest(*p)
        p_plot_mean = []
        p_plot_std = []
        for i in itertools.izip_longest(*p):
            i = np.asarray(i)
            i = i[i != np.array(None)]
            p_plot_mean.append( np.mean(i))
            p_plot_std.append( np.std(i))
        p_plot_mean = np.asarray(p_plot_mean)
        p_plot_std = np.asarray(p_plot_std)
        ax.plot(np.log10(range(1, len(p_plot_mean )+1)),p_plot_mean, label='Gens in seed bank = ' + str(int(1 / (c/M))), lw = 2, color = colors[j])
        ax.fill_between(np.log10(range(1, len(p_plot_mean )+1)), p_plot_mean+p_plot_std, p_plot_mean-p_plot_std, facecolor=colors[j], alpha=0.5)

        if max(g) > maxg:
            maxg = max(g)
    ax.set_xlim([0, np.log10(maxg)])
    ax.set_ylim([0, 1])
    ax.legend(loc='upper left')
    ax.set_xlabel('Generations, ' + r'$log_{10}$', fontsize=20)
    ax.set_ylabel('Frequency of favored allele', fontsize=20)
    fig.savefig(mydir + 'figs/SweepDorm.png', \
        bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

#sweepDormFig()
