import datetime

import pandas as pd
import numpy as np

import init
import params
from Bicluster import bicluster

import multiprocessing as mp
## see: http://matthewrocklin.com/blog/work/2013/12/05/Parallelism-and-Serialization/
## and https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-pythons-multiprocessing-pool-ma
## and https://stackoverflow.com/questions/19984152/what-can-multiprocessing-and-dill-do-together
##import pathos.multiprocessing as mp

## Note these par funcs need to be in globals.py (and hence global) to have global access to all
##    data -- I don't know how to change this right now.

def do_something_par( items, func, threads=4 ):
    pool = mp.Pool(processes=threads)              # start 4 worker processes
    ## Need to send, e.g. a tuple (1, counts_g) if fill_all_scores_par() took multiple args
    out = pool.map(func, items)
    pool.terminate()
    return out

def fill_all_cluster_scores_par( clusters, threads=4 ):
    clusters = do_something_par( clusters, fill_all_scores_par, threads=threads )
    clusters = {clusters[i].k: clusters[i] for i in xrange(len(clusters))}  # convert back to map
    return clusters

## Keep everything global so it doesn't need to be sent to each child.
def fill_all_scores_par( k ):
    global clusters, all_genes, ratios, string_net, counts_g
    clust = clusters[k]
    ##if np.all(np.invert(clust.changed)):
    ##    return clust
    clust.fill_all_scores(all_genes, ratios, string_net, counts_g, ratios.columns.values)
    return clust

iter = 1

ratios, genome_seqs, anno, op_table, string_net, allSeqs_fname, all_bgFreqs, all_genes = \
    init.pynkey_init(params.organism, params.k_clust, params.ratios_file)

# Save all pynkey code for safe keeping
pynkey_code = {}
try:
    pynkey_code = init.load_pynkey_code()
except:
    print 'Must have "import"ed this module rather than run from main.py.'
    
startTime = datetime.datetime.now()
print str(startTime)

clusters = init.init_biclusters( ratios, params.k_clust, 'kmeans+random' );
        
counts_g = bicluster.get_all_cluster_row_counts( clusters, all_genes )

print 'DONE WITH INITIALIZATION!'
endTime = datetime.datetime.now()
print str(endTime)
print str(endTime - startTime) + ' seconds since initialization'
        
stats_df = pd.DataFrame()

