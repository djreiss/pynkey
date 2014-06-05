## MAIN PROGRAM
import datetime

import pandas as pd

import init
import params
from Bicluster import bicluster
import funcs
import utils as ut
import scores
import floc

def run_pynkey(iter):
    global organism, n_iters, clusters, all_genes, ratios, string_net, stats_df

#     n_no_improvements = 0
#     for i=iter:n_iters
#         iter = i
#         (clusters, n_improvements, stats_tmp) = 
    floc.do_floc( clusters, iter, all_genes, ratios, string_net )
#         println( @sprintf( "%.3f", (time() - startTime)/60 ), " minutes since initialization" )
#         stats_df = rbind( stats_df, stats_tmp )
#         write_table( "output/$(organism)_stats.tsv", stats_df )
#         if isfile( "DO_SAVE" )  ## save cluster info for temporary examination of clusters (via Rscripts/clusters.R)
#             warn( "Writing out clusters to output/$(organism)_clusters.tsv" )
#             clusters_tab = clusters_to_dataFrame(clusters);
#             write_table("output/$(organism)_clusters.tsv", clusters_tab)
#         end
#         if n_improvements <= 0 n_no_improvements += 1 else n_no_improvements = 0; end
#         if iter > n_iters/2 && n_no_improvements > 5 break; end
#     end
    return iter
# end

if __name__ == '__main__':
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

    ##clusters = bicluster.fill_all_cluster_scores_par( clusters, all_genes, ratios, string_net, 
    ##                                 ratios.columns.values, 4 )
    clusters = bicluster.fill_all_cluster_scores( clusters, all_genes, ratios, string_net, ratios.columns.values )

    print 'DONE WITH INITIALIZATION!'
    endTime = datetime.datetime.now()
    print str(endTime)
    print str(endTime - startTime) + ' seconds since initialization'
        
    stats_df = pd.DataFrame()

    ##iter = run_pynkey(iter) ## Note this function can be run like this to restart from current iter

    print 'DONE!'
    endTime = datetime.datetime.now()
    print str(endTime)
    print str(endTime - startTime) + ' seconds since initialization'

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
