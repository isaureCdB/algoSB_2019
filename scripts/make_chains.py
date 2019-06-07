#!/usr/bin/env python3

"""
    make_chains.py
    Assembles all chains based on a connectivity graph, and a maximum meanrank
    In addition, for every chain, the exact overlap RMSD is computed from the preatom and postatom coordinates,
        (and the l-RMSD toward the bound form if L-RMSD for each fragment is provided)

    Copyright 2015-2017 Sjoerd de Vries, Isaure Chauvot de Beauchene, TUM
"""

import sys, json, argparse, json
import numpy as np
from math import log

########################
parser =argparse.ArgumentParser(description=__doc__,
formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument('graph', help="connectivity graph (in JSON format), \
                            as generated by connect.py or meanrank_filter.py")
parser.add_argument('--meanrank', help="maximum (geometric) meanrank of the chains")
parser.add_argument('--preatoms', nargs='+', help="preatom coordinates for each fragment,\
                                                in numpy (.npy) format")
parser.add_argument('--postatoms', nargs='+', help="postatom coordinates for each fragment,\
                                                in numpy (.npy) format")
parser.add_argument('--lrmsd', nargs='+', help="lrmsd of each fragment (optional)")
parser.add_argument('--maxchains', help="Maximal number of chains", default=100000)

args = parser.parse_args()
########################

def write_chain(indices, sum_overlap_msds, meanrank):
    '''  write the chain in output file'''
    global count
    meanrank = "%.3f" % meanrank
    o = np.sqrt(sum_overlap_msds / (len(indices) - 1) )
    o = "%.3f" % o
    lr = []
    if args.lrmsd is not None:
        for inr, i in enumerate(indices):
            lr.append(lrmsds[inr][i])
        lrms = np.sqrt(sum([v*v for v in lr])/len(lr))
        lrms = "%.3f" % lrms
        print("#indices", lrms, meanrank, o, end=' ')
    else:
        print("#indices", meanrank, o, end=' ')
    for i in indices: print(i, end=' ')
    for l in lr: print(l, end=' ')
    print('')
    count += 1
    if not count%10000: print(count, file=sys.stderr)
    if count >= maxchains:
        sys.exit("Finding more than %i chains.Use a lowest meanrank or a highest maxchains"%maxchains)

def walk(pos, curr, indices, sum_overlap_msds, curr_meanrank):
    '''
    Do the chain building, iteratively.
      pos = position of fragment in the chain
      curr = pose index in the pool for the considered fragment
      indices = which pose taken for each of the upstream fragments
    '''
    ind = clusters[pos][curr]
    new_indices = indices + (ind,)
    curr_meanrank += log(ind)
    # if the chain already has a too high score:
    if curr_meanrank > meanrank_threshold: return
    overlap_msd = 0 # msd of the overlapping parts in the chain
                    # (how well the poses fit spatially together)
    if pos > 0:
        pre = preatoms[pos-1][indices[-1]-1]
        post = postatoms[pos][ind-1]
        d = post - pre
        overlap_msd = (d*d).sum()/pre.shape[0]
    new_sum_overlap_msds = sum_overlap_msds + overlap_msd
    if pos == nfrag - 1: # if we reached the end of the chain
        meanrank = np.exp(curr_meanrank/nfrag)
        write_chain(new_indices, new_sum_overlap_msds, meanrank)
        return
    for target in interactions[pos][curr]:
        # take the next connected pose to the current pose
        walk(pos+1, target, new_indices, new_sum_overlap_msds, curr_meanrank)


if __name__ == "__main__":
    # The json file contains the connectivity graph
    # This is the output either of connect.py,
    # or of meanrank_filter.py applied to the output of connect.py.
    # It contains the list of connected poses/custers
    # for each pair of consecutive fragments
    tree = json.load(open(args.graph))

    # Nb of fragments in the chain
    nfrag = tree["nfrags"]
    assert nfrag == len(args.preatoms), (nfrag, len(args.preatoms))
    assert nfrag == len(args.postatoms), (nfrag, len(args.postatoms))
    if args.lrmsd is not None:
        assert nfrag == len(args.lrmsd), (nfrag, len(args.lrmsd))

    # User-defined maximal score/energy for the whole chain
    max_meanrank = float(args.meanrank)
    meanrank_threshold = log(max_meanrank) * nfrag

    maxchains = int(args.maxchains)

    # np.array of poses coordinates
    postatoms = [np.load(f) for f in args.postatoms]
    preatoms = [np.load(f) for f in args.preatoms]
    for arr in (preatoms, postatoms):
        for anr, a in enumerate(arr):
            ncoor = int(a.shape[1] / 3)
            arr[anr] = a.reshape(a.shape[0], ncoor, 3)

    lrmsds = []
    if args.lrmsd is not None:
            # extract the ligand-rmsd of the poses computed by ATTRACT
            # This is not used inthe chain-building, only written in
            # the output file
        for f in args.lrmsd:
            lrmsd = {}
            lnr = 0
            for l in open(f):
                lnr += 1
                ll = l.split()
                if len(ll) != 2: continue
                k = ll[0]
                if k == "l-RMSD": k = lnr
                else: k = int(k)
                value = float(ll[1])
                lrmsd[k] = value
            lrmsds.append(lrmsd)

    # clusters = list-of-lists-of-Cluster objects, 1 list per frag
    # extracted from the tree
    clusters = []
    for cnr, tclus in enumerate(tree["clusters"]):
        clus = []
        for tcc in tclus:
        # tclus is a dictionary
        # for now, limit ourselves to fully collapsed trees:
        # there is either only one pose per cluster or all poses
        # of the cluster are within 0.1 A (deredundant cutoff, see connect.py)
            assert tcc["radius"] == 0 # clustering radius = 0 (collapsed tree)
            rank = np.array(tcc["ranks"], dtype="int")[0]
            clus.append(rank)
        clusters.append(clus)

    # for each frag, interactions = list of the connected downstream poses
    # extracted from the tree
    interactions = [{} for l in range(nfrag-1)]
    for cnr, tinter in enumerate(tree["interactions"]):
        inter = interactions[cnr]
        for source, target in tinter:
            #(source, target) is a pair of connected poses.
            # Store the connections with the next fragment.
            # Count number of connections with previous fragment.
            if source not in inter: inter[source] = []
            inter[source].append(target)

    count = 0 # count number of chains

    print("#header <mean (root-mean-sq) ligand rmsd> <mean (geometric mean) rank>  <rms-overlap-rmsd> <indices> <ligand rmsds>")
    f = clusters[0]
    for ff in range(len(f)):
        walk(0, ff, (), 0, 0)

    print("found %i chains"%count, file=sys.stderr)