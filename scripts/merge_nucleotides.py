#!/usr/bin/env python3

'''
average atom coordinates on overlapping parts of nucl steps
usage: chain2rna.py <chains_file> --nat <nat1 nat2 nat3>
        --coor <step1.npy step2.npy step3.npy>
'''

import sys, argparse, numpy as np
from npy import npy2to3, npy3to2

#######################################
parser =argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('chains_file', help="txt file, one line per structure, \
                    each line containing the list of fragment indices (from 1)\
                     at each position of the structure, with space separator.")
parser.add_argument('--nat', help="list of the numbers of atoms per nucleotide", nargs='+', type=int)
# avoid loading twice the same motif.npy, which is time consuming
parser.add_argument('--coor', help="coordinates of one step for all structures,\
                    (one file per step), in npy format with shape (nstruct, ncoor)\
                    with ncoor = 3*nat", nargs='+')
parser.add_argument('--outp')
args = parser.parse_args()
#######################################

chainfile = args.chains_file # txt file. One chain per line, one index per column
try:
    steps = np.loadtxt(chainfile, dtype=int) - 1 ### changed 9-02_17:30
    if len(steps.shape) == 1:
        steps = np.reshape(steps,(1,steps.shape[0]))
except:
    print("no np.loadtxt")
    cc = [ [int(i)-1 for i in l.split()] for l in open(chainfile)]
    steps = np.array(cc)

nfrag = steps.shape[1]
nat = args.nat
print(nat, file=sys.stderr)

assert nfrag == len(nat) + 1, (nfrag, len(nat))

coor = [ npy2to3(np.load(i)) for i in args.coor ] # one np.array per step

outp = args.outp
assert outp != args.chains_file, "ERROR: output is same as input"

#initialise merged structure
max_atom = sum([n.shape[1] for n in coor])  #max nb of atoms in final chain
rna = np.zeros( (len(steps), max_atom, 3) )
count = 0

#First atoms unchanged
len_step = coor[0].shape[1]   #nb of atoms in 1st step
n = len_step - nat[0]         #nb of atom to not merge in 1st step
rna[:,:n] = coor[0][steps[:, 0], :n]
count += n

for i in range(nfrag-1):
    #Merge overlapping atoms
    print("merge nucl %i"%(i+2), file=sys.stderr)
    coor1 = coor[i][steps[:, i], -nat[i]: ]      # last atoms of previous step
    coor2 = coor[i+1][steps[:, i+1], :nat[i] ]   # first atoms of next step
    rna[:, count:count+nat[i]] = 0.5*(coor1+coor2)
    count += nat[i]
    if i < nfrag-2:
        #add non-overlapping atoms of next step
        len_step = coor[i+1].shape[1]
        n = len_step - nat[i] - nat[i+1]
        print(coor[i+1].shape)
        rna[:, count:count+n] = coor[i+1][steps[:, i+1], nat[i]:nat[i]+n]
        count += n

#Last atoms unchanged
len_step = coor[-1].shape[1]   #nb of atoms in 1st step
n = len_step - nat[-1]         #nb of atom to not merge in 1st step
rna[:,count:count+n] = coor[-1][steps[:, -1], -n:]
count += n

rna = rna[:, :count]

np.save(outp, rna)
