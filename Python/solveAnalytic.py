from __future__ import division
import numpy as np
import pandas as pd
import os
from scipy import integrate
import  matplotlib.pyplot as plt

mydir = os.path.expanduser('~/github/DormSweep/')

# Definition of parameters
s = 0.1
r = 0.7
d = 0.3
def dP_dt(P, t=0):
    ''' Return the change in frequency of the favored allele for active and dormant
    individuals. '''
    '''P is an array with the allele frequencies of active (P[0]) and dormant (P[0]) individuals'''
    return np.array( [(s * P[0] * (1-P[0]))  - (d * P[0]) + (r * P[1]) ,
                        (-P[1] * r ) + (d * P[0])])


P_f0 = np.array([ 0. ,  0.])
P_f1 = np.array([ 1, d/r ])
all(dP_dt(P_f0) == np.zeros(2) ) and all(dP_dt(P_f1) == np.zeros(2)) # => True

def d2P_dt2(P, t=0):
    """ Return the Jacobian matrix evaluated in X. """
    return np.array([[ -d + s - (2 * P[0]) * s , r  ],
                  [ d,  -r ] ])

#print d2P_dt2(P_f0)
#print d2P_dt2(P_f1)

def plotDynamics():
    fig = plt.figure()
    t = np.linspace(0, 1000,  1000)              # time
    P0 = np.array([0.01, 0])                     # initials conditions: 10 rabbits and 5 foxes
    P, infodict = integrate.odeint(dP_dt, P0, t, full_output=True)
    infodict['message']                     # >>> 'Integration successful.'
    active, dormant = P.T
    plt.plot(t, active, 'r-', label='Active')
    plt.plot(t, dormant  , 'b-', label='Dormant')
    plt.grid()
    plt.legend(loc='best')
    plt.xlabel('Time [generations]', fontsize=20)
    plt.ylabel('Allele frequency', fontsize=20)
    plt.title('Selective sweep of a favored allele')
    plt.savefig(mydir + 'figs/testDyn.png')


plotDynamics()
