from __future__ import division
import pandas as pd
import os, argparse, math, itertools
import numpy as np

mydir = os.path.expanduser('~/github/DormSweep/')

def cleanData(N = 1000, M = 1000, s = 0.1, r = 100):
    #Cs = [1, 10, 100]
    #Cs = [100]
    maxg = 0
    pops = ['N', 'M']
    Cs = np.logspace(0, 4, num = 100, base=10.0)
    for pop in pops:
        df = pd.DataFrame()
        for j, c in enumerate(Cs):
            df_j = pd.read_csv(mydir + 'data/sweeps_dorm/' + pop + '_sweep_N_' + str(N) + '_M_' + str(M) + '_c_' + \
                str(c) + '_s_' + str(s) + '_r_' + str(r) + '.txt', header = None)
            times = df_j.count(axis=1)
            times.columns = str(c)
            df[str(c)] = times

        df.to_csv(mydir +'data/T_fix/T_fix_' +  pop + '_sweep_N_' + str(N) + '_M_' + str(M) + \
            '_s_' + str(s) + '_r_' + str(r) + '.txt' , sep = ' ', index = False)


cleanData()
