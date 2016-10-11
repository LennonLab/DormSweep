from __future__ import division
import numpy as np
import pandas as pd
import os

mydir = os.path.expanduser("~/github/DormSweep/")

def sweepSim(N, s, reps = 1000):
    count = 0
    freqs_list = []
    while count < reps:
    # B = selected allele
        g = 0
        B_count = 1
        b_count = N-1
        freqs = [1/N]
        gs = [1]
        g = 1
        while B_count < N and B_count != 0:
            b_b =  (1 - (B_count/N)) * (1 - (B_count/N))
            b_B = (B_count/N) * (1 - (B_count/N))
            B_b = (B_count/N) * (1 - (B_count/N)) * (1-s)
            B_B = ((B_count/N) * (B_count/N)) +  ( (B_count/N) * (1 - (B_count/N)) * s)
            pop = list(np.random.choice(2, N, p=[b_b + B_b, b_B + B_B]))
            B_count = pop.count(1)
            b_count = pop.count(0)
            freqs.append(B_count / N)
            g += 1
            gs.append(g)
            if B_count == 0:
                break
        if B_count != N:
            continue
        else:
            print N, s, count
            freqs_list.append(freqs)
            count +=1
    #print freqs_list
    df = pd.DataFrame(freqs_list)
    df.to_csv(mydir + 'data/sweep_' + str(s) + '.txt', header=False, index = False)


def multipleS(N, reps = 100):
    Ss = [0.001, 0.01, 0.1]
    for s in Ss:
        sweepSim(N, s = s, reps = reps)

multipleS(1000)
