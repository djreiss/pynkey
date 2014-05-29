## MAIN PROGRAM
import datetime

import params
import init

iter = 1

ratios, genome_seqs, anno, op_table, string_net, allSeqs_fname, all_bgFreqs, all_genes = \
    init.pynkey_init(params.organism, params.k_clust, params.ratios_file)

# Save all pynkey code for safe keeping
pynkey_code = init.load_pynkey_code()

startTime = datetime.datetime.now()
print str(startTime)

clusters = init.init_biclusters( ratios, params.k_clust, 'kmeans+random' );
# if nprocs() > 1 clusters = fill_all_cluster_scores_parallel( clusters, true, true );
# else clusters = fill_all_cluster_scores( clusters, true, true ); end
# println( @sprintf( "%.3f", (time() - startTime)/60 ), " minutes since initialization" )

# stats_df = DataFrame()
# run_junkey() ## Note this function can be run like this to restart from current iter

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
