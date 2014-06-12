## MAIN PROGRAM
import datetime

import pandas as pd
from params import n_iters,nthreads

def run_pynkey(iter):
    n_no_improvements = 0
    for i in range(iter,n_iters):
        iter = i
        globals.clusters, n_improvements, stats_tmp = floc.run( globals.clusters, iter, globals.all_genes, \
                                                            globals.ratios, globals.string_net, globals.startTime )
        ##println( @sprintf( "%.3f", (time() - startTime)/60 ), " minutes since initialization" )
        globals.stats_df = globals.stats_df.append( stats_tmp )
        ##write_table( "output/$(organism)_stats.tsv", stats_df )
        ##if isfile( "DO_SAVE" ): ## save cluster info for temp. examination of clusters
        ##    warn( "Writing out clusters to output/$(organism)_clusters.tsv" )
        ##    clusters_tab = clusters_to_dataFrame(clusters);
        ##    write_table("output/$(organism)_clusters.tsv", clusters_tab)
        if n_improvements <= 0:
            n_no_improvements += 1
        else:
            n_no_improvements = 0
        if iter > n_iters/2 and n_no_improvements > 5:
            break
    globals.iter = iter + 1
    return globals.iter

if __name__ == '__main__':
    import globals ## this does the initialization
    from Bicluster import fill_all_cluster_scores_par
    import floc

    #clusters = fill_all_cluster_scores( clusters, all_genes, ratios, string_net, ratios.columns.values )    
    ## weird - if I move this to globals.py, then it gets locked up.
    globals.clusters = fill_all_cluster_scores_par(globals.clusters, threads=nthreads)

    # NOTE: run_pynkey() which calls floc.get_floc_scores_all() fills all the cluster scores at the beginning    
    globals.iter = run_pynkey(globals.iter) ## Note this function can be run like this to restart from current iter

    print 'DONE!'
    globals.endTime = datetime.datetime.now()
    print str(globals.endTime)
    print str(globals.endTime - globals.startTime) + ' seconds since initialization'

# kInd = 1;
# if organism == "Hpy" || organism == "Eco" kInd = clusters_w_func("flagell", clusters)[1]; end
# b = clusters[kInd]
# #b = re_meme_bicluster( b );
# print(join(b.meme_out, '\n'))

# try mkdir( "output" ); end
# ## Right now this fails on the clusters if we use gzopen (open works fine, though)
# save_jld( "output/$(organism)_out.jldz", (organism, k_clust, ratios, genome_seqs, anno, op_table, string_net, 
#                                    allSeqs_fname, all_bgFreqs, startTime, endTime,
#                                    all_genes, iter, n_iters, distance_search, distance_scan, motif_width_range,
#                                    clusters, stats_df, junkey_code) )

# clusters_tab = clusters_to_dataFrame(clusters);
# write_table("output/$(organism)_clusters.tsv", clusters_tab) ## for examination of clusters in R (via Rscripts/clusters.R)

# tmp = get_cluster_row_counts(clusters);
# println( sum(tmp.==0), " genes in no clusters" )
# println( sum(tmp.==maximum(tmp)), " genes in ", maximum(tmp), " clusters" )

# println( @sprintf( "%.3f", (endTime - startTime)/60 ), " minutes since initialization" )


# #genes = rownames(ratios)[clusters[kInd].rows] ##rows]
# #seqs = get_sequences(genes);
# #@time gibbs_out = gibbs_site_sampler(seqs[:,2])     ## run gibbs sampler on most "flagellar-enriched" cluster
# #@time gibbs_out2 = gibbs_site_sampler(seqs, gibbs_out["pssm"])

# end; ## myid() == 1
