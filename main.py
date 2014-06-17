## MAIN PROGRAM
import datetime
import os.path
import warnings

import pandas as pd
import numpy as np
from params import n_iters,nthreads,organism
import funcs

def run_pynkey(iter):
    n_no_improvements = 0
    for i in range(iter,n_iters):
        iter = i
        glb.iter = iter
        clusters, n_improvements, stats_tmp = floc.run( glb.clusters, iter, glb.all_genes, glb.ratios, glb.string_net )
        glb.clusters = clusters
        ##println( @sprintf( "%.3f", (time() - startTime)/60 ), " minutes since initialization" )
        glb.stats_df = glb.stats_df.append( stats_tmp )
        if os.path.exists( 'DO_SAVE' ): ## save stats, and cluster info for temp. examination of clusters
            warnings.warn( 'Writing out clusters to output/%s_clusters.tsv' % organism )
            glb.stats_df.to_csv( 'output/%s_stats.tsv' % organism, sep='\t', na_rep='NA' )
            clusters_tab = clusters_to_dataFrame(clusters)
            clusters_tab.to_csv( 'output/%s_clusters.tsv' % organism, sep='\t', na_rep='NA' )

        n_no_improvements = n_no_improvements+1 if n_improvements <= 0 else 0
        n_changed = np.nansum( [np.sum(clust.changed) for clust in glb.clusters.values()] )
        if n_changed <= 10 and iter > n_iters/2 and n_no_improvements > 5:
            break
    return glb.iter

def finish():
    from Bicluster import bicluster

    print 'DONE!'
    glb.endTime = datetime.datetime.now()
    print str(glb.endTime)
    print str(glb.endTime - glb.startTime) + ' seconds since initialization'

    kInd = 1
    if organism == 'Hpy' or organism == 'Eco':
        kInd = funcs.clusters_w_func('flagell', glb.clusters, glb.anno)[0]

    b = glb.clusters[kInd]
    print b.meme_out

## Right now this fails on the clusters if we use gzopen (open works fine, though)
# save_jld( "output/$(organism)_out.jldz", (organism, k_clust, ratios, genome_seqs, anno, op_table, string_net, 
#                                    allSeqs_fname, all_bgFreqs, startTime, endTime,
#                                    all_genes, iter, n_iters, distance_search, distance_scan, motif_width_range,
#                                    clusters, stats_df, junkey_code) )

    clusters_tab = funcs.clusters_to_dataFrame(glb.clusters)
    clusters_tab.to_csv( 'output/%s_clusters.tsv' % params.organism, sep='\t', na_rep='NA' )

    tmp = np.array( bicluster.get_all_cluster_row_counts( glb.clusters, glb.all_genes ).values() )
    print np.sum(tmp==0), 'genes in no clusters'
    print np.sum(tmp==np.max(tmp)), 'genes in', np.max(tmp), 'clusters'

# println( @sprintf( "%.3f", (endTime - startTime)/60 ), " minutes since initialization" )

# #genes = rownames(ratios)[clusters[kInd].rows] ##rows]
# #seqs = get_sequences(genes);
# #@time gibbs_out = gibbs_site_sampler(seqs[:,2])     ## run gibbs sampler on most "flagellar-enriched" cluster
# #@time gibbs_out2 = gibbs_site_sampler(seqs, gibbs_out["pssm"])


if __name__ == '__main__':
    import globals as glb ## this does the initialization
    from Bicluster import fill_all_cluster_scores_par
    import floc

    try:
        os.mkdir( 'output' )
    except:
        print 'Cannot mkdir ./output/'

    #clusters = fill_all_cluster_scores( clusters, all_genes, ratios, string_net, ratios.columns.values )    
    ## weird - if I move this to glb.py, then it gets locked up.
    glb.clusters = fill_all_cluster_scores_par(glb.clusters, threads=nthreads)

    # NOTE: run_pynkey() which calls floc.get_floc_scores_all() fills all the cluster scores at the beginning    
    glb.iter = run_pynkey(glb.iter) ## Note this function can be run like this to restart from current iter

    finish()
