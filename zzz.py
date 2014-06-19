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

### Plotting scores to figure out whats going on
from matplotlib import pyplot as plt
import utils

counts_g = bicluster.get_all_cluster_row_counts( clusters, all_genes )
counts_g2 = np.array( [ val for k,val in counts_g.items() ] )
score_g = bicluster.get_cluster_row_count_scores( counts_g )
score_g2 = np.array( [ score_g[k] for k in counts_g ] )

utils.setup_text_plots(fontsize=8, usetex=True)
plt.scatter( counts_g2, score_g2 )
