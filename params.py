import sys
import numpy as np

print 'importing params'

nthreads = 1 ## 1 to not parallelize; None to automatically use all available processors

#def undefined(var):
#    var not in vars() and var not in globals()

#if undefined('organism'):
organism = 'Hpy'
#organism = 'Eco'
#organism = 'Sce'
#organism = 'Mpn'

print organism

#if undefined('k_clust'):
if organism == 'Hpy':
    k_clust = 75
elif organism == 'Eco':
    k_clust = 450
elif organism == 'Sce':
    k_clust = 600
elif organism == 'Mpn':
    k_clust = 75

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

max_network_weight = 0.9
max_motif_weight =   1.8

avg_genes_per_cluster = 22
avg_clusters_per_gene = 1.3

## allow more updates if there are more clusters??? Tuned to k_clust/2 for Hpy (where k_clust is 75) --
##     may need additional tuning; e.g. for eco (k_clust=450), k_clust/2 is too high
max_improvements_per_iter = round(min(k_clust/2, 100))

distance_search = np.array(distance_search)
distance_scan = np.array(distance_scan)
motif_width_range = np.array(motif_width_range)
