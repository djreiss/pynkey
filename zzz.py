##%run main.py

### reload run data/clusters
import init
init.init( 'output/Hpy.pkl' )

from globals import *
import funcs
import floc

if False:
    scores_all = floc.get_scores_all(clusters, iter, all_genes, ratios, string_net)
    scores_all2 = floc.get_scores_best(scores_all)
    #%timeit -n1 ord = floc.rnd_bubblesort( scores_all2['combined'].values.copy() ) ##, scores_all2.shape[0]*2 )
    #%timeit -n1 ord = floc.rnd_bubblesort2( scores_all2['combined'].values.copy() ) ##, scores_all2.shape[0]*2 )
    #%timeit -n1 ord = floc.rnd_bubblesort3( scores_all2['combined'].values.copy() ) ##, scores_all2.shape[0]*2 )


### testing the parallelized biclusterings
from globals import *
from params import *
import Bicluster
tmp = Bicluster.re_meme_all_clusters_par(clusters)

### testing weaved residue funcs
import weaved

b = clusters[0]
print funcs.matrix_residue( ratios.ix[b.rows, b.cols] )
rats = ratios.ix[b.rows, b.cols].values
print weaved.fast_resid(rats)


%timeit funcs.matrix_residue( rats, weaveIt=True )
## 10000 loops, best of 3: 35 us per loop
%timeit funcs.matrix_residue( rats, weaveIt=False )
## 1000 loops, best of 3: 489 us per loop

## NOTE that the lookups are now taking up 90% of the time!
%timeit funcs.matrix_residue( ratios.ix[b.rows, b.cols], weaveIt=True )
## 1000 loops, best of 3: 687 us per loop
%timeit funcs.matrix_residue( ratios.ix[b.rows, b.cols], weaveIt=False )
## 1000 loops, best of 3: 1.1 ms per loop
%timeit ratios.ix[b.rows, b.cols].values
## 1000 loops, best of 3: 600 us per loop


### How to debug the Bicluster.py code:
reload(Bicluster)
from Bicluster import bicluster
bb=Bicluster.bicluster(b.k,b.rows,b.cols)


for bb in clusters.values():
    rats = ratios.ix[ bb.rows, bb.cols ]
    rats2 = ratios.ix[ bb.rows, ratios.columns.values[~np.in1d(np.array(ratios.columns.values,str), bb.cols)] ]
    print np.nanmean( np.abs( rats.values ) ), np.nanmean( np.abs( rats2.values ) )


### Plotting scores to figure out whats going on
from matplotlib import pyplot as plt
import utils

counts_g = bicluster.get_all_cluster_row_counts( clusters, all_genes )
counts_g2 = np.array( [ val for k,val in counts_g.items() ] )
score_g = bicluster.get_cluster_row_count_scores( counts_g )
score_g2 = np.array( [ score_g[k] for k in counts_g ] )

utils.setup_text_plots(fontsize=8, usetex=True)
plt.scatter( counts_g2, score_g2 )



## using networkx
import networkx as nx
net=nx.read_weighted_edgelist('Eco/STRING_511145.tsv.gz')
## get sum of all weights of edges between given nodes! 
rows = glb.clusters[0].rows
## this uses networkx and takes 1.92 ms:
np.sum(nx.get_edge_attributes(net.subgraph(rows),'weight').values()) 
## this uses pandas and takes 226 ms:
np.sum(glb.string_net[ glb.string_net[[0,1]].isin(rows).all(1) ].weight)/2.0
## even this takes longer (27.7 ms + 4.33 ms) = 32.03 ms:
net2 = glb.string_net.ix[rows]
np.sum( net2[ net2[[0,1]].isin(rows).all(1) ].weight)/2.0
## Using networkx to do this on all 4000 genes is still (27.7+4.33*4000)/(1.92*4000) or about 2.25x faster

## how about this? no -- 114 ms -- why? and does it account for weights?
nx.average_node_connectivity(net.subgraph(rows))


## using igraph? can't read in gzipped file, so uncompress it then:
import igraph as ig
G=ig.Graph.Read_Ncol('Eco/STRING_511145.tsv',weights=True)
## throws an error if any rows are not in the network
## This is 4.05 ms
np.sum(G.induced_subgraph(rows[np.in1d(rows,G.vs['name'])]).es['weight'])/2.0
r=rows[np.in1d(rows,G.vs['name'])]
## This is 250 us -- so if we could pre-filter all rows into only those that are in the network, this is fastest
np.sum(G.induced_subgraph(r).es['weight'])/2.0

