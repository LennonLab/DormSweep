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
    df = pd.DataFrame(freqs_list)
    df.to_csv(mydir + 'data/sweep_' + str(s) + '.txt', header=False, index = False)


def sweepSimDorm(N, M, c, s, reps = 1000):
    count = 0
    N_freqs_list = []
    M_freqs_list = []
    while count < reps:
    # B = selected allele
        Mb = M
        MB = 0
        NB = 1
        Nb = N-1
        N_freqs = [(NB/N)]
        M_freqs = [(MB/M)]
        gs = [1]
        g = 1
        K = N/M
        d = (c* K) / N
        r = (c  ) / M
        while (NB + MB) < (N + M) and (NB + MB) != 0:
            # transistion probs https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2206092/
            Nb_Nb =  (1-  (NB/N)) * (1-  (NB/N)) * (1-d)
            Nb_NB = (NB/N) * (1- (NB/N)) * (1-d)
            NB_Nb = (NB/N) * (1- (NB/N)) * (1-d) * (1-s)
            NB_NB = ( (NB/N) * (NB/N) * (1-d) ) + ((NB/N) * (1- (NB/N)) * s * (1-d) )

            Nb_Mb = (1 - (NB/N)) * d
            Nb_MB = 0
            NB_MB = (NB/N) * d
            NB_Mb = 0

            Mb_Nb = (1 - (MB/M)) * r
            Mb_NB = 0
            MB_NB = (MB/M) * r
            MB_Nb = 0

            Mb_Mb = (1- (MB/M)) * (1-r)
            Mb_MB = 0
            MB_Mb = 0
            MB_MB = (MB/M) * (1-r)

            # chose b and B
            #print Nb_Nb + NB_Nb + Mb_Nb + MB_Nb + Nb_NB + NB_NB + Mb_NB + MB_NB
            #print Nb_Mb + NB_Mb + Mb_Mb + MB_Mb + Nb_MB + NB_MB + Mb_MB + MB_MB
            pop_N = list(np.random.choice(2, N, p=[Nb_Nb + NB_Nb + Mb_Nb + MB_Nb, \
                Nb_NB + NB_NB + Mb_NB + MB_NB]))
            pop_M = list(np.random.choice(2, M, p=[Nb_Mb + NB_Mb + Mb_Mb + MB_Mb, \
                Nb_MB + NB_MB + Mb_MB + MB_MB]))

            Nb = pop_N.count(0)
            NB = pop_N.count(1)
            Mb = pop_M.count(0)
            MB = pop_M.count(1)
            N_freqs.append((NB / N) )
            M_freqs.append((MB / M) )
            g += 1
            gs.append(g)
            if (NB + MB) == 0:
                break
        if (NB + MB)  != (N + M):
            continue
        else:
            print N, M, c, s, count
            N_freqs_list.append(N_freqs)
            M_freqs_list.append(M_freqs)
            count +=1
    N_df = pd.DataFrame(N_freqs_list)
    M_df = pd.DataFrame(M_freqs_list)
    N_df.to_csv(mydir + 'data/N_sweep_N_' + str(N) + '_M_' + str(M) + '_c_' + \
        str(c) + '_s_' + str(s) + '_r_' + str(reps) + '.txt', header=False, index = False)
    M_df.to_csv(mydir + 'data/M_sweep_N_' + str(N) + '_M_' + str(M) + '_c_' + \
        str(c) + '_s_' + str(s) + '_r_' + str(reps) + '.txt', header=False, index = False)

def multipleS(N, reps = 100):
    Ss = [0.001, 0.01, 0.1]
    for s in Ss:
        sweepSim(N, s = s, reps = reps)

def multipleSDorm(N, M, s = 0.1, reps = 100):
    #Ss = [0.001, 0.01, 0.1]
    #Cs = [1, 10, 100]
    #Cs = np.logspace(0, 4, num = 100, base=10.0)
    Cs = [10000]
    for c in Cs:
        print c
        sweepSimDorm(N = N, M = M, c = c,  s = s, reps = reps)

multipleSDorm(N = 1000, M = 10000)
