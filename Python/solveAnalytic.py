from __future__ import division
import numpy as np
import pandas as pd
import os, itertools
from scipy import integrate
import  matplotlib.pyplot as plt

mydir = os.path.expanduser('~/github/DormSweep/')


# Jacobian class
class modelMatrix():
    def __init__(self, s, r, d):
        self.r = r
        self.s = s
        self.d = d

    def dP_dt(self, P, t=0):
        ''' Return the change in frequency of the favored allele for active and dormant
        individuals. '''
        '''P is an array with the allele frequencies of active (P[0]) and dormant (P[0]) individuals'''
        return np.array( [(self.s * P[0] * (1-P[0]))  - (self.d * P[0]) + (self.r * P[1]) ,
                            (-P[1] * self.r ) + (self.d * P[0])])

    def d2P_dt2(P, t=0):
        """ Return the Jacobian matrix evaluated in X. """
        return np.array([[ -self.d + self.s - (2 * P[0]) * self.s , self.r  ],
                      [ self.d,  -self.r ] ])

def plotDynamics(N = 1000):
    fig = plt.figure()
    s = 0.1
    t = np.linspace(0, 10000,  10000)              # time
    P0 = np.array([1/N, 0])
    rs_and_ds = [0.001, 0.01, 0.1]
    colors = ['#87CEEB', '#FFA500', '#FF6347']
    for i, r_and_d in enumerate(rs_and_ds):
        P, infodict = integrate.odeint(modelMatrix(s,r_and_d,r_and_d).dP_dt, \
            P0, t, full_output=True)
        #infodict['message']                     # >>> 'Integration successful.'
        active, dormant = P.T
        pop_mean = (active + dormant) / 2
        plt.plot(np.log10(t), pop_mean, lw = 2, label= r'$m_{N,M}, m_{M,N} = $' + str(r_and_d), color = colors[i])
    plt.grid()
    plt.ylim(0,1)
    plt.xlim(0,4)
    plt.legend(loc='upper left', fontsize = 12)
    plt.xlabel('Time (generations), ' + r'$log_{10}$', fontsize=20)
    plt.ylabel('Frequency of favored allele', fontsize=20)
    #plt.title('Selective sweep of a favored allele')
    plt.savefig(mydir + 'figs/testDyn.png', bbox_inches = "tight", \
        pad_inches = 0.4, dpi = 600)

def plotNs(N = 1000):
    fig = plt.figure()
    s = 0.1
    t = np.linspace(0, 100000,  100000)
    #N = np.logspace(2, 10000,  10000)
    P0 = np.array([1/N, 0])
    Cs = np.logspace(0, 4, num = 100, base=10.0)
    T_fix_A = []
    T_fix_D = []
    T_fix_mean = []
    for C in Cs:
        P, infodict = integrate.odeint(modelMatrix(s,1/C,1/C).dP_dt, \
            P0, t, full_output=True)
        #infodict['message']                     # >>> 'Integration successful.'
        active, dormant = P.T
        #active = np.ndarray.tolist(active)
        #dormant = np.ndarray.tolist(dormant)
        #print active[active >= float(1)]
        T_fix_A_m = np.argmax(active>=float(0.997))
        T_fix_D_m = np.argmax(dormant>=float(0.997))
        T_fix_A.append(T_fix_A_m)
        T_fix_D.append(T_fix_D_m)
        T_fix_mean.append((T_fix_A_m + T_fix_D_m) /2 )

    pops = ['N', 'M']
    fig, ax = plt.subplots()
    colors = ['#87CEEB',  '#FF6347']
    pop_type = {'N': 'Active', 'M': 'Dormant'}
    for i, pop in enumerate(pops):
        if pop == 'N':
            y = T_fix_A
        else:
            y = T_fix_D
        plt.plot(np.log10(Cs), y, lw = 2, label= pop_type[pop], color = colors[i])

    plt.plot(np.log10(Cs), T_fix_mean, lw = 2, color = 'g')

    ax.set_xlabel('Average time in seed bank, ' + r'$log_{10}$', fontsize=20)
    plt.grid()
    plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, \
           ncol=2, mode="expand", borderaxespad=0.)
    plt.axvline(x = np.log10(1000), linewidth=2, color='darkgrey',ls='--')
    ax.set_ylim([0,35000])

    ax.set_ylabel( r'$T_{fix}$', fontsize=20)
    fig.savefig(mydir + 'figs/Tfix_N_M_model.png', \
        bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

def plot_sim_v_model(N=1000):
    fig, ax = plt.subplots()
    colors = ['#87CEEB', '#FFA500', '#FF6347']
    s = 0.1
    t = np.linspace(0, 10000,  10000)
    cs = [1, 10, 100]
    P0 = np.array([1/N, 0])
    for x, c in enumerate(cs):
        P, infodict = integrate.odeint(modelMatrix(s,1/c,1/c).dP_dt, \
            P0, t, full_output=True)
        active, dormant = P.T
        pop_mean = (active + dormant) / 2
        sim = pd.read_csv(mydir + 'data/sweeps_dorm/sweep_N_1000_M_1000_c_' + \
            str(c) + '_s_0.1.txt', header = None)
        g = []
        p = []
        sim = sim.fillna(1.0)
        for index, row in sim.iterrows():
            g.extend(row[~row.isnull()].index.values)
            p.append(row[~row.isnull()].values)

        p_plot_mean = []
        p_plot_std = []
        for i in itertools.izip_longest(*p):
            i = np.asarray(i)
            i = i[i != np.array(None)]
            p_plot_mean.append( np.mean(i))
            p_plot_std.append( np.std(i))
        p_plot_mean = np.asarray(p_plot_mean)
        p_plot_std = np.asarray(p_plot_std)
        print c, len(pop_mean), len(p_plot_std)
        pop_mean = pop_mean[:len(p_plot_mean)]
        print c, len(pop_mean), len(p_plot_std)
        ax.plot(pop_mean, p_plot_mean, lw = 2, color = colors[x], \
            label='Time in seed bank = ' + str(int(1 / (c/1000))))
        ax.fill_between(pop_mean, p_plot_mean+p_plot_std, p_plot_mean-p_plot_std, facecolor=colors[x], alpha=0.5)
    plt.grid()
    plt.plot([0, 1],[0, 1], linewidth=2, color='darkgrey',ls='--')
    ax.set_xlim([0,1])
    ax.set_ylim([0,1])
    ax.set_xlabel('Allele frequency from logistic model', fontsize=20)
    ax.set_ylabel('Allele frequency from simulation', fontsize=20)
    ax.legend(loc='upper left', fontsize = 12)
    fig.savefig(mydir + 'figs/sim_vs_model.png', \
        bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()



def plot_sim_v_model_active(N=1000):
    fig, ax = plt.subplots()
    s = 0.1
    t = np.linspace(0, 100000,  100000)
    c=10.235310219
    P0 = np.array([1/N, 0])
    P, infodict = integrate.odeint(modelMatrix(s,1/c,1/c).dP_dt, \
        P0, t, full_output=True)
    active, dormant = P.T
    N_sim = pd.read_csv(mydir + 'data/sweeps_mac/N_sweep_N_1000_M_10000_c_10.235310219_s_0.1.txt', header = None)
    M_sim = pd.read_csv(mydir + 'data/sweeps_mac/M_sweep_N_1000_M_10000_c_10.235310219_s_0.1.txt', header = None)
    g_N = []
    p_N = []
    g_M = []
    p_M = []

    N_sim = N_sim.fillna(1.0)
    for index, row in N_sim.iterrows():
        g_N.extend(row[~row.isnull()].index.values)
        p_N.append(row[~row.isnull()].values)

    p_N_plot_mean = []
    p_N_plot_std = []
    for i in itertools.izip_longest(*p_N):
        i = np.asarray(i)
        i = i[i != np.array(None)]
        p_N_plot_mean.append( np.mean(i))
        p_N_plot_std.append( np.std(i))
    p_N_plot_mean = np.asarray(p_N_plot_mean)
    p_N_plot_std = np.asarray(p_N_plot_std)
    active = active[:len(p_N_plot_mean)]
    ax.plot(active, p_N_plot_mean, lw = 2, color = '#87CEEB')
    ax.fill_between(active, p_N_plot_mean+p_N_plot_std, p_N_plot_mean-p_N_plot_std, facecolor='#87CEEB', alpha=0.5)
    plt.grid()
    fig.savefig(mydir + 'figs/sim_vs_model_active.png', \
        bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.close()

plotDynamics()
#plot_sim_v_model()
#plot_sim_v_model()
