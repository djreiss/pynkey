import sys

import numpy as np

print 'importing params'

organism = 'Hpy'
#organism = 'Eco'
#organism = 'Sce'
#organism = 'Mpn'

n_iters = 100
k_clust = 0

## these are defaults (for small microbes like Halo):
distance_search = [-150,+10]
distance_scan =   [-250,+20]
motif_width_range = [6,24]

max_resid_weight =   1.0
max_motif_weight =   0.5 ##2.0 ##38.0 ##1.8
max_network_weight = 0.1 ##0.5 ##1.0 ##29.0  ##0.9
max_volume_weight =  0.1 ## 1.0 ## 3.0
max_clusters_per_gene_weight = 0.01 ##1.0 ##3.0 ##0.1
max_column_weight = 0.2

avg_genes_per_cluster = 22
avg_clusters_per_gene = 1.8 ##1.3

nthreads = None ## 1 to not parallelize; None to automatically use all available processors

## set random seed here!
import time
import random
random_seed = int(time.clock()*1e6) ## 1403806278241227 ## this is good one for hpy
print 'RANDOM SEED:', random_seed
np.random.seed( random_seed )
random.seed( random_seed ) ## I dont know if I use non-numpy random, but do it here just in case

ratios_file = None
max_improvements_per_iter = None

## if 'all' then use all gene names possible (vs. only those listed in the expression data)
all_genes_option = 'only_in_expression_data' ## 'all' 

output_dir = './output/'

def parse_args():
    global organism, k_clust, nthreads
    global max_resid_weight, max_motif_weight, max_network_weight, max_volume_weight, max_clusters_per_gene_weight
    options = None

    from optparse import OptionParser
    parser = OptionParser(usage='%prog [options]')
    parser.add_option('-o', '--organism', dest='organism', 
                        help='input organism [default: %default]', metavar='ORGANISM', default='Hpy')
    parser.add_option('--out', dest='output_dir', default=output_dir,
                        help='output directory [default: %default]', metavar='OUT')
    parser.add_option('-k', '--k', dest='k', type=int, default=0, ## default is different for given organism
                        help='number of clusters [default: %default]', metavar='KCLUST')
    parser.add_option('-r', '--r0', '--resid_weight', dest='resid_weight', type=float, default=max_resid_weight, 
                        help='maximum residual weight [default: %default]', metavar='WEIGHT')
    parser.add_option('-m', '--m0', '--mot_weight', dest='mot_weight', type=float, default=max_motif_weight, 
                        help='maximum motif weight [default: %default]', metavar='WEIGHT')
    parser.add_option('-n', '--n0', '--net_weight', dest='net_weight', type=float, default=max_network_weight,
                        help='maximum network weight [default: %default]', metavar='WEIGHT')
    parser.add_option('-v', '--v0', '--volume_weight', dest='volume_weight', type=float, default=max_volume_weight,
                        help='maximum volume weight [default: %default]', metavar='WEIGHT')
    parser.add_option('-g', '--g0', '--clusters_per_gene_weight', dest='clusters_per_gene_weight', type=float, default=max_clusters_per_gene_weight,
                        help='maximum clusters_per_gene weight [default: %default]', metavar='WEIGHT')
    parser.add_option('--nthreads', dest='nthreads', type=int, default=8,
                        help='maximum number of threads used [default: %default]', metavar='N')
    parser.add_option('--pylab') ## ignore this if run from ipython
    
    options, args = parser.parse_args()
    ##return(options,args)

    try: ## fails in ipython if called via 'import params'
        options,args = parse_args()
        print options.organism
        if options is None:   ## --help called (throws error, above)
            sys.exit(0)
    except:
        None

    if options is not None:

        try:
            nthreads = options.nthreads
        except:
            nthreads = None

        try:
            organism = options.organism
        except:
            None ##organism = 'Eco'

        try:
            output = options.output_dir
        except:
            None
            
        try:
            k_clust = options.k
        except:
            None

        try:
            max_resid_weight = options.resid_weight
        except:
            None

        try:
            max_motif_weight = options.mot_weight
        except:
            None

        try:
            max_network_weight = options.net_weight
        except:
            None

        try:
            max_volume_weight = options.volume_weight
        except:
            None

        try:
            max_clusters_per_gene_weight = options.clusters_per_gene_weight
        except:
            None

def init_args():
    global organism, k_clust, nthreads, ratios_file, output_dir
    global max_resid_weight, max_motif_weight, max_network_weight, max_volume_weight, max_clusters_per_gene_weight
    global motif_width_range, distance_search, distance_scan, max_improvements_per_iter

    print 'ORGANISM:', organism

    if k_clust is None or k_clust == 0:
        if organism == 'Hpy' or organism == 'Mpn':
            k_clust = 75
        elif organism == 'Eco':
            k_clust = 450
        elif organism == 'Sce':
            k_clust = 600

    print k_clust

    try:
        os.mkdir( output_dir )
    except:
        print 'Cannot mkdir ' + output_dir

    #if undefined('ratios'):
    ratios_file = './' + organism + '/ratios.tsv.gz'
    print ratios_file

    if organism == 'Eco':
        motif_width_range = [6,30]
    elif organism == 'Sce':
        distance_search = [-250,+30]
        distance_scan =   [-500,+50]
        motif_width_range = [6,18]

    print distance_search
    print distance_scan
    print motif_width_range

    ## allow more updates if there are more clusters??? Tuned to k_clust/2 for Hpy (where k_clust is 75) --
    ##     may need additional tuning; e.g. for eco (k_clust=450), k_clust/2 is too high
    max_improvements_per_iter = round(min(k_clust/2, 100))

    distance_search = np.array(distance_search)
    distance_scan = np.array(distance_scan)
    motif_width_range = np.array(motif_width_range)

