from optparse import OptionParser

import numpy as np

print 'importing params'

options = None

def parse_args():
    parser = OptionParser()
    parser.add_option('-o', '--organism', dest='organism',
                        help='input organism [default: %default]', metavar='ORGANISM', default='Hpy')
    parser.add_option('-k', '--k', dest='k', type=int, default=75,
                        help='number of clusters [default: %default]', metavar='KCLUST')
    parser.add_option('-r', '--r0', '--resid_weight', dest='resid_weight', type=float, default=1.0, 
                        help='maximum residual weight [default: %default]', metavar='RESID_WEIGHT')
    parser.add_option('-m', '--m0', '--mot_weight', dest='mot_weight', type=float, default=2.0, 
                        help='maximum motif weight [default: %default]', metavar='MOT_WEIGHT')
    parser.add_option('-n', '--n0', '--net_weight', dest='net_weight', type=float, default=0.5,
                        help='maximum network weight [default: %default]', metavar='NET_WEIGHT')
    parser.add_option('--nthreads', dest='nthreads', type=int, default=8,
                        help='maximum number of threads used [default: %default]', metavar='NTHREADS')
    #parser.add_option('--pylab') ## ignore this if run from ipython
    
    options, args = parser.parse_args()
    return(options,args)

try: ## fails in ipython if called via 'import params'
    options,args = parse_args()
    print options.organism
except:
    None

nthreads = None ## 1 to not parallelize; None to automatically use all available processors
try:
    nthreads = options.nthreads
except:
    nthreads = None

## set random seed here!
import time
import random
random_seed = int(time.time()*1e6)
print 'RANDOM SEED:', random_seed
np.random.seed( random_seed )
random.seed( random_seed ) ## I dont know if I use non-numpy random, but do it here just in case

try:
    organism = options.organism
except:
    #organism = 'Hpy'
    organism = 'Eco'
    #organism = 'Sce'
    #organism = 'Mpn'

print 'ORGANISM:', organism

try:
    k_clust = options.k
except:
    if organism == 'Hpy' or organism == 'Mpn':
        k_clust = 75
    elif organism == 'Eco':
        k_clust = 450
    elif organism == 'Sce':
        k_clust = 600

#if undefined('ratios'):
ratios_file = './' + organism + '/ratios.tsv.gz'

print ratios_file

n_iters = 100

## these are defaults (for small microbes like Halo):
distance_search = [-150,+10]
distance_scan =   [-250,+20]
motif_width_range = [6,24]

if organism == 'Eco':
    motif_width_range = [6,30]
elif organism == 'Sce':
    distance_search = [-250,+30]
    distance_scan =   [-500,+50]
    motif_width_range = [6,18]

print distance_search
print distance_scan
print motif_width_range

max_resid_weight =   1.0
max_motif_weight =   2.0 ##38.0 ##1.8
max_network_weight = 0.5 ##1.0 ##29.0  ##0.9
max_volume_weight =  0.3 ## 3.0
max_clusters_per_gene_weight = 1.0 ##3.0 ##0.1
max_column_weight = 0.2

try:
    max_resid_weight = options.resid_weight
except:
    max_resid_weight = 1.0

try:
    max_motif_weight = options.mot_weight
except:
    max_motif_weight = 1.0 # 2.0

try:
    max_network_weight = options.net_weight
except:
    max_network_weight = 0.25 # 0.5

avg_genes_per_cluster = 22
avg_clusters_per_gene = 1.3

## allow more updates if there are more clusters??? Tuned to k_clust/2 for Hpy (where k_clust is 75) --
##     may need additional tuning; e.g. for eco (k_clust=450), k_clust/2 is too high
max_improvements_per_iter = round(min(k_clust/2, 100))

distance_search = np.array(distance_search)
distance_scan = np.array(distance_scan)
motif_width_range = np.array(motif_width_range)

## if 'all' then use all gene names possible (vs. only those listed in the expression data)
all_genes_option = 'only_in_expression_data' ## 'all' 

