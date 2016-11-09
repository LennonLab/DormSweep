from __future__ import division
import numpy as np
import pandas as pd
import os
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



def plotDynamics():
    fig = plt.figure()
    s = 0.1
    t = np.linspace(0, 10000,  10000)              # time
    P0 = np.array([0.01, 0])
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


plotDynamics()
