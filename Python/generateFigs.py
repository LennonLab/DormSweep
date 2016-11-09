from __future__ import division
import pandas as pd
import  matplotlib.pyplot as plt
import os, argparse, math, itertools
import numpy as np
import scipy as sp
import scipy.stats

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
    ax.set_ylabel('Frequency of allele', fontsize=20)
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
        df = pd.read_csv(mydir + 'data/sweeps_dorm/sweep_N_' + str(N) + '_M_' + str(M) + '_c_' + \
            str(c) + '_s_' + str(s) + '.txt', header = None)
        df = df.fillna(1.0)
        #print df

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
        ax.plot(np.log10(range(1, len(p_plot_mean )+1)),p_plot_mean, label='Time in seed bank = ' + str(int(1 / (c/M))), lw = 2, color = colors[j])
        ax.fill_between(np.log10(range(1, len(p_plot_mean )+1)), p_plot_mean+p_plot_std, p_plot_mean-p_plot_std, facecolor=colors[j], alpha=0.5)

        if max(g) > maxg:
            maxg = max(g)
    ax.set_xlim([0, np.log10(maxg)])
    ax.set_ylim([0, 1])
    plt.grid()
    ax.legend(loc='upper left', fontsize = 12)
    ax.set_xlabel('Time (generations), ' + r'$log_{10}$', fontsize=20)
    ax.set_ylabel('Frequency of favored allele', fontsize=20)
    fig.savefig(mydir + 'figs/SweepDorm.png', \
        bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def mean_confidence_interval(data, confidence=0.95):
    a = 1.0*np.array(data)
    n = len(a)
    m, se = np.mean(a), scipy.stats.sem(a)
    h = se * sp.stats.t._ppf((1+confidence)/2., n-1)
    return m, m-h, m+h

def TfixFig(M = 10000):
    pops = ['N', 'M']
    fig, ax = plt.subplots()
    colors = ['#87CEEB',  '#FF6347']
    pop_type = {'N': 'Active', 'M': 'Dormant'}
    IN_N = pd.read_csv(mydir + 'data/T_fix/T_fix_N_sweep_N_1000_M_10000_s_0.1_r_100.txt', sep = ' ')
    IN_M = pd.read_csv(mydir + 'data/T_fix/T_fix_M_sweep_N_1000_M_10000_s_0.1_r_100.txt', sep = ' ')
    df_add = IN_N.add(IN_M, fill_value=0)
    df = df_add.divide(2, axis=0)
    Cs = df.columns.values.astype(float)
    timeInSb = 1 / (Cs / M)
    means = []
    std = []

    for column in df:
        data = df[column].values.astype(float)
        means.append(np.mean(data))
        std.append(np.std(data))
    #means = np.log10(np.asarray(means))
    means = np.asarray(means)
    std = np.asarray(std)
    plt.axvline(x = np.log10(1000), linewidth=2, color='darkgrey',ls='--')

    ax.plot(np.log10(timeInSb), means,  lw = 2, color = '#87CEEB')
    plt.grid()
    ax.fill_between(np.log10(timeInSb), means+std, means-std, facecolor='#87CEEB', alpha=0.5)

    ax.set_xlabel('Average time in seed bank, ' + r'$log_{10}$', fontsize=20)
    ax.set_ylabel( r'$T_{fix}$', fontsize=20)
    fig.savefig(mydir + 'figs/Tfix.png', \
    bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



def TfixFigActDorm(M = 10000):
    pops = ['N', 'M']
    fig, ax = plt.subplots()
    colors = ['#87CEEB',  '#FF6347']
    pop_type = {'N': 'Active', 'M': 'Dormant'}
    for i, pop in enumerate(pops):
        IN = pd.read_csv(mydir + 'data/T_fix/T_fix_' + str(pop) + '_sweep_N_1000_M_10000_s_0.1_r_100.txt', sep = ' ')
        Cs = IN.columns.values.astype(float)

        timeInSb = 1 / (Cs / M)
        means = []
        m_plus_hs = []
        m_minus_hs = []
        std = []
        for column in IN:
            data = IN[column].values.astype(float)
            #m, m_plus_h, m_minus_h = mean_confidence_interval(data, confidence=0.95)
            means.append(np.mean(data))
            std.append(np.std(data))
            #m_plus_hs.append(m_plus_h)
            #m_minus_hs.append(m_minus_hs)
            #print ms, m_plus_hs, m_minus_hs
        means = np.asarray(means)
        std = np.asarray(std)
        #print Cs
        #print ms
        #print ms
        #print len(ms), len(timeInSb)
        #p_plot_mean = np.asarray(ms)
        #p_plot_mean_plus_h = np.asarray(m_plus_hs)
        #p_plot_mean_minus_h = np.asarray(m_minus_hs)
        #print len(p_plot_mean), len(timeInSb)

        ax.plot(np.log10(timeInSb), means, label= pop_type[pop], lw = 2, color = colors[i])
        #ax.fill_between(np.log10(range(1, len(timeInSb )+1)), means, m_minus_hs, facecolor=colors[i], alpha=0.5)
        ax.fill_between(np.log10(timeInSb), means+std, means-std, facecolor=colors[i], alpha=0.5)

    ax.set_xlabel('Average time in seed bank, ' + r'$log_{10}$', fontsize=20)
    plt.grid()
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, \
           ncol=2, mode="expand", borderaxespad=0.)
    plt.axvline(x = np.log10(1000), linewidth=2, color='darkgrey',ls='--')

    ax.set_ylabel( r'$T_{fix}$', fontsize=20)
    fig.savefig(mydir + 'figs/Tfix_N_M.png', \
        bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

TfixFigActDorm()
sweepDormFig()
TfixFig()
