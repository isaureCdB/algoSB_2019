#!/usr/bin/env python3

'''
average atom coordinates on overlapping parts of fragments in chains
usage: chain2rna.py <chains> <nfrag> <nat> <motifs> <motif[1-m].npy> <rna.npy>
'''

import sys, argparse, numpy as np
from npy import npy2to3, npy3to2

#######################################
parser =argparse.ArgumentParser(description=__doc__,
                        formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('chains_file', help="txt file, one line per structure, \
    each line containing the list of fragment indices (from 1) at each position\
    of the structure, with space separator. \n \
    Output from make_chains.py or meanrank_filter.py. ex: UUU-5frag-2A.chains")
parser.add_argument('--nfrag', help="Nb of fragments per chain", type=int)
parser.add_argument('--nat', help="Nb of atoms per nucleotide", nargs='+', type=int)
parser.add_argument('--motifs', help="index of the motif in motifs list: \
                    ex: 1 1 2 for AUU + UUU + UUU", nargs='+', type=int)
# avoid loading twice the same motif.npy, which is time consuming
parser.add_argument('--npy', help="npy of ach motif, once per motif: \
                    ex: AUU.npy UUU.npy ", nargs='+')
parser.add_argument('--outp')
args = parser.parse_args()
#######################################

chainfile = args.chains_file # txt file. 1 chain per line, 1 pose index per column
try:
    chains = np.loadtxt(chainfile, dtype=int) - 1 ### changed 9-02_17:30
    if len(chains.shape) == 1:
        chains = np.reshape(chains,(1,chains.shape[0]))
except:
    print("no np.loadtxt")
    cc = [ l for l in open(chainfile).readlines() if not l.startswith("#header") ]
    c = cc[0].split()
    if c[1][0] == "#indices":
        c = [k[1:] for k in c]
    for i in range(len(c)):
        try:
            int(c[i])
            break
        except:
            continue
    nfrag = int((len(c)-3)*0.5)
    chains = np.array([ [int(i)-1 for i in l.split()[3:3+nfrag]] for l in cc ] )

if args.nfrag is None:
    nfrag = len(chains[0])
else:
    nfrag = args.nfrag
    assert nfrag == len(chains[0])

nat = args.nat[:nfrag+2] # 7666 for AUUU
m = [int(i)-1 for i in args.motifs[:nfrag]] # motifs: 1 1 2 for AUU+UUU+UUU
#assert nfrag == len(nat) - 2
assert len(m) == nfrag, (m, nfrag)

struct = [ npy3to2(np.load(i)) for i in args.npy ] # one np.array per motif
outp = args.outp
assert outp != args.chains_file, "outp is same as input"

c = [ 3*n_atoms for n_atoms in nat ]
for i in range(nfrag):
    assert struct[m[i]].shape[1] == sum(c[i:i+3]), (i, struct[m[i]].shape , sum(c[i:i+3]))

rna = np.zeros( (len(chains),3*sum(nat)) )

#First nucleotide unchanged
rna[:,:c[0]] = struct[m[0]][ chains[:,0] ,:c[0]]

#Merge pairs of overlapping atoms in 1st overlap
coor1 = struct[m[0]][chains[:,0], c[0]:sum(c[:2])]
coor2 = struct[m[1]][chains[:,1], :c[1]]
rna[:, c[0]:c[0]+c[1]]  = 0.5*(coor1 + coor2)

#Merge triplets of overlapping atoms
for i in range(2, nfrag):
    coor1 = struct[m[i-2]][chains[:,i-2], sum(c[i-2:i]):]
    coor2 = struct[m[i-1]][chains[:,i-1], c[i-1]:-c[i+1] ]
    coor3 = struct[m[i-0]][chains[:,i-0], :c[i]]
    rna[:, sum(c[:i]):sum(c[:i+1])] = (coor1+coor2+coor3)/3

#Merge pairs of overlapping atoms in last overlap
coor1 = struct[m[-3]][chains[:,-2], -c[-2]:]
coor2 = struct[m[-1]][chains[:,-1], c[-3]:sum(c[-3:-1])]
rna[:, sum(c[:-2]):sum(c[:-1])] = 0.5*( coor1 + coor2 )

#Last nucleotide unchanged
rna[:, sum(c[:-1]):sum(c)] = struct[m[-1]][ chains[:,-1] ,-c[nfrag-1]:]

np.save(outp, rna)
