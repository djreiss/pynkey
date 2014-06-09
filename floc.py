# Use FLOC algorithm to update clusters

import pandas as pd
import numpy as np
from numpy import random as rnd
from numpy import nan as NA

import utils as ut
import scores
from Bicluster import bicluster,fill_all_cluster_scores_par

print 'importing floc'

## Get gain scores for all possible row/col moves
## Default: allow best 5 row- and 20 col- moves per bicluster instead of all of them!
def get_floc_scores_all(clusters, iter, all_genes, ratios, string_net, max_row=9999, max_col=9999):  ## 5, 20)
    ## First, collate move scores into a single DataFrame for stochastic sorting
    ## Use pd.DataFrame for scores rather than matrix

    ##bicluster.fill_all_cluster_scores(clusters, all_genes, ratios, string_net, ratios.columns.values)
    clusters = fill_all_cluster_scores_par(clusters, threads=15)

    counts_g = bicluster.get_all_cluster_row_counts( clusters, all_genes ) ## Clusters per gene counts - precompute
    tmp = [ clusters[k].get_floc_scoresDF_rows(iter, all_genes, ratios, counts_g) for k in clusters ]
    scoresDF_r = pd.concat(tmp) ## This is much faster than mapreduce() above!!!
    del tmp

    ## STANDARDIZE ROWS AND COLUMNS SEPARATELY!!!   -- is this really what we want?
    tmp1 = scores.get_combined_scores( scoresDF_r, iter, ratios )
    tmp1 = ut.sdize_vector( tmp1 )
    scoresDF_r['combined'] = pd.Series(tmp1)

    tmp = [ clusters[k].get_floc_scoresDF_cols(iter, ratios) for k in clusters ]
    scoresDF_c = pd.concat(tmp) ## This is much faster than mapreduce()
    del tmp

    weight_r, weight_n, weight_m, weight_c, weight_v, weight_g = scores.get_score_weights( iter, ratios )

    tmp1 = scores.get_combined_scores( scoresDF_c, iter, ratios )
    tmp1 = ut.sdize_vector( tmp1 ) * weight_c
    scoresDF_c['combined'] = pd.Series(tmp1)

    if max_row < 9999:
        tmp_r = scoresDF_r.groupby( 'k' )
        tmp_r = dict(list(tmp_r))
        shrunk_r = [pd.DataFrame() for i in range(len(tmp_r))]
        for i in range(len(tmp_r)):
            shrunk_r[i] = tmp_r[i].sort( 'combined' )[ 0:max_row ]
        scoresDF_r = pd.concat( shrunk_r )

    if max_col < 9999:
        tmp_c = scoresDF_c.groupby( 'k' )
        tmp_c = dict(list(tmp_c))
        shrunk_c = [pd.DataFrame() for i in range(len(tmp_c))]
        for i in range(len(tmp_c)):
            shrunk_c[i] = tmp_c[i].sort( 'combined' )[ 0:max_col ]
        scoresDF_c = pd.concat( shrunk_c )

    all_scores = pd.concat([scoresDF_r, scoresDF_c])
    return all_scores

# ## Instead of scores for ALL possible moves, make a matrix of n_best best scores for each row/col
def get_floc_scores_best( all_scores, n_best_row=3, n_best_col=3 ):
    df_r = all_scores[ all_scores.is_row_col == 'r' ]
    if n_best_row < 9999:
        tmp = df_r.groupby( 'row_col' ) ## Cool!
        dfs_r = {i:pd.DataFrame() for i,df in tmp}
        for r,tmp2 in tmp:
            tmp2a = tmp2.sort( 'combined' )[ 0:n_best_row ]
            dfs_r[r] = tmp2a

    df_c = all_scores[ all_scores.is_row_col == 'c' ]
    if n_best_col < 9999:
        tmp = df_c.groupby( 'row_col' ) ## Cool!
        dfs_c = {i:pd.DataFrame() for i,df in tmp}
        for r,tmp2 in tmp:
            tmp2a = tmp2.sort( 'combined' )[ 0:n_best_col ]
            dfs_c[r] = tmp2a

    scores_out = pd.concat( [ pd.concat( dfs_r ), pd.concat( dfs_c ) ] )
    return scores_out

## Native code: 119s
def rnd_bubblesort( scores, Nrepeats=None ):
    lsc = len(scores)
    if Nrepeats == None: ## is None:
        Nrepeats = lsc * 2
    ord = np.arange(lsc)
    rnd.shuffle(ord) ## start w/ random order
    tmp = ut.minmax(scores)
    R = 2.0 * ( tmp[1]-tmp[0] ) ## Denominator of value to compute prob. from
    the_max = tmp[1]
    n = lsc - 1
    sc = scores.copy()
    sc[ np.isnan(sc) ] = the_max ## replace NaN with maximum score
    n_switches = 0
    for i in xrange(Nrepeats): ## TBD: this double loop can be sped up with weave!!!
        o1 = o2 = g1 = g2 = p = 0.
        rnds = rnd.rand(n)
        for j in xrange(n):
            o1 = ord[j]
            o2 = ord[j+1]
            g1 = sc[o1]
            g2 = scores[o2]
            if g1 == g2 and g2 == the_max:
                continue
            p = 0.5 + ( g1 - g2 ) / R ## compute prob of switching
            if rnds[j] < p: ##rnd.rand() < p: ## switch???
                ord[j] = o2
                ord[j+1] = o1
                n_switches += 1
        if i % 1000 == 1:
            print( i, n_switches, Nrepeats )
    return ord

## This is jit-c'd via LVVM using numba - takes about 127 secs
## seems that it is not really jit-ed
from numba import float64, int64, jit
rnd_bubblesort2 = jit( int64[:]( double[:], int64 ) )( rnd_bubblesort )

## this is weaved (using c++) -- takes about 3.04 seconds
import scipy.weave
def rnd_bubblesort3( scores, Nrepeats=None ):
    lsc = len(scores)
    if Nrepeats is None:
        Nrepeats = lsc * 2
    ord = np.arange(lsc)
    rnd.shuffle(ord) ## start w/ random order
    tmp = ut.minmax(scores)
    R = 2.0 * ( tmp[1]-tmp[0] ) ## Denominator of value to compute prob. from
    the_max = tmp[1]
    n = lsc - 1
    sc = scores.copy()
    sc[ np.isnan(sc) ] = the_max ## replace NaN with maximum score
    switchesN = np.array([0]) ## count the number of switches. Not really necessary
    for i in xrange(Nrepeats):
        rnds = rnd.rand(n)
        code = \
"""
  for ( int j=0; j<n; j++ ) {
      int o1 = ord(j), o2 = ord(j+1);
      double g1 = sc(o1), g2 = sc(o2);
      if ( g1 == g2 && g2 == (double)the_max ) continue;
      double p = 0.5 + ( g1 - g2 ) / (double)R; // compute prob of switching
      if ( rnds(j) < p ) { // switch???
          ord(j) = o2;
          ord(j+1) = o1;
          switchesN(0) ++;
      }
    }
 """
        scipy.weave.inline(code,['ord','n','sc','R','switchesN','rnds','the_max'],
                           type_converters=scipy.weave.converters.blitz)
        if i % 1000 == 1:
            print i, switchesN[0], Nrepeats
    return ord

# ## TODO: add max_improvements param (to prevent really fast optimization at beginning before motifing turns on)
def floc_update(clusters, iter, all_genes, ratios, string_net, max_no_improvements=25):

    scores_all = get_floc_scores_all(clusters, iter, all_genes, ratios, string_net)

    scores_all2 = get_floc_scores_best(scores_all)

    ## For FLOC, identify the BEST score for each row/column, then bubble sort those to get the order in
    ##    which to perform the moves.
    ## Note this is wrong right now - it sorts ALL k scores for each row/col. 
    ##  Need to just use the BEST score for each row/col and then bubblesort these.
    ord = rnd_bubblesort2( scores_all2['combined'].values ) ##, n_sort_iter) 
    print scores_all2.ix[ord,:].head(); print scores_all2.ix[ord,:].tail()

#     new_clusters = saved_clusters = copy_clusters( clusters, true, false ); ## make a copy for updating
#     (weight_r, weight_n, weight_m, weight_c, weight_v, weight_g) = get_score_weights() ## DONE: don't need to update n or m scores if their weights are 0

#     all_resids = float32( [ clusters[k].resid for k=1:k_clust ] )
#     all_dens_string = float32( [ clusters[k].dens_string for k=1:k_clust ] )
#     all_meanp_meme = float32( [ clusters[k].meanp_meme for k=1:k_clust ] )
#     all_volumes = zeros( Float32, length(all_resids) )
#     mean_resid = nanmean(all_resids)
#     mean_dens_string = nanmean(all_dens_string) ## Should usually increase, sometimes get worse
#     mean_meanp_meme = nanmean(all_meanp_meme) ## Should usually decrease, sometimes get worse
#     counts_gene = get_cluster_row_counts( clusters ) ## Clusters per gene score
#     #scores_gene = get_cluster_row_count_scores( counts_gene )
#     #mean_score_gene = nanmean(scores_gene)
#     no_improvements = 0 ## Number of updates in which we haven't seen an improvement
#     n_improvements = 0; n_tries = 0
#     for i=ord ## Update bicusters -- should "store" the one with the best mean resid during the updates
#         sc = scores2[i,:]
#         kUpd = sc["k"][1]
#         cc::bicluster = copy_cluster(new_clusters[kUpd], true, true)
#         row_col = sc["row_col_ind"][1]
#         if sc["is_row_col"][1] == 'r'
#             ##if sc["is_in"][1] && length(cc.rows) <= min_rows continue; end ## Don't remove if we're already at the min. Now: use volume score instead
#             cc.rows = sc["is_in"][1] ? remove(cc.rows, row_col) : [cc.rows, row_col]
#             cc.changed[1] = true
#         elseif sc["is_row_col"][1] == 'c'
#             cc.cols = sc["is_in"][1] ? remove(cc.cols, row_col) : [cc.cols, row_col]
#             cc.changed[2] = true
#         end
#         cc.resid = bicluster_residue( cc, ratios )
#         all_resids[kUpd] = cc.resid
#         if sc["is_row_col"][1] == 'r'  ## Only update network/motif scores if it's a row, duh!
#             if abs(weight_n) > 0            ## only compute if need to
#                 cc.dens_string = bicluster_network_density( cc, string_net )
#                 all_dens_string[kUpd] = cc.dens_string
#             end
#             if weight_m > 0            ## only compute if need to
#                 cc.meanp_meme = bicluster_meme_pval( cc )
#                 all_meanp_meme[kUpd] = cc.meanp_meme
#             end
#             ## TODO: incorporate score_g (nclust per gene) and score_v (ngene per clust) into combined scores
# #             if weight_g > 0 
# #                 counts_gene[ row_col ] += (sc["is_in"][1] ? -1 : 1) ## Clusters per gene score
# #                 mean_score_gene = nanmean(get_cluster_row_count_scores(cc, counts_gene)) ## Clusters per gene score
# #                 counts_gene[ row_col ] -= (sc["is_in"][1] ? -1 : 1) ## revert it back
# #             end
#         end
#         ##all_scores[kUpd] = get_combined_score( cc.resid, cc.dens_string, cc.meanp_meme, 0.0f0 )
#         ##all_scores = get_combined_scores( sdize_vector( all_resids ), sdize_vector( all_dens_string ), 
#         ##                                 sdize_vector( all_meanp_meme ), all_volumes )
#         ## Need to normalize these deltas vs. the stddev of the corresponding scores in the "score" dataframe!
#         score_delta = get_combined_score( (nanmean(all_resids) - mean_resid) / nansd(all_resids), 
#                                          (nanmean(all_dens_string) - mean_dens_string) / nansd(all_dens_string),
#                                          (nanmean(all_meanp_meme) - mean_meanp_meme) / nansd(all_meanp_meme), 
#                                          0.0f0, 0.0f0 ) ##(nanmean(scores_gene) - mean_score_gene) / nansd(scores_gene) )
#         new_clusters[kUpd] = cc
#         n_tries += 1
#         ##println(nanmean(all_scores)," ",mean_scores)
#         if score_delta <= 0 ## mean bicluster scores improved, so store the entire cluster stack
#             n_improvements += 1
#             mean_resid = nanmean(all_resids) ## Should usually decrease, sometimes get worse
#             mean_dens_string = nanmean(all_dens_string) ## Should usually increase, sometimes get worse
#             mean_meanp_meme = nanmean(all_meanp_meme) ## Should usually decrease, sometimes get worse
#             ##mean_scores = nanmean(all_scores) ## Should usually decrease, sometimes get worse
#             ##saved_clusters[kUpd] = cc ##copy_cluster(cc, true)
#             saved_clusters = copy_clusters( new_clusters, false, false ) ## make a copy to keep the best update
#             output = @sprintf( "%d %.4f %.4f %.4f %.4f %d %c %s %d %.4f %.4f %.4f %d %d %d",
#                               i, mean_resid, mean_dens_string, mean_meanp_meme, score_delta, ##mean_scores,
#                               Base.int(kUpd), sc["is_row_col"][1], sc["is_in"][1]?"remove":"add", row_col, 
#                               cc.resid, clusters[kUpd].resid, cc.resid-clusters[kUpd].resid, 
#                               n_tries, n_improvements, no_improvements)
#             println( output )
#             no_improvements = 0
#         else
#             no_improvements += 1
#         end
#         if no_improvements > max_no_improvements break; end
#     end
#     (saved_clusters, n_improvements, n_tries, scores2[ord,:])
# end

def do_floc(clusters, iter, all_genes, ratios, string_net):
#     ##clusters = fill_cluster_scores(clusters) ## dont need this since each cluster's scores are updated in floc_update
#     ## allow more updates if there are more clusters??? Tuned to k_clust/2 for Hpy (where k_clust is 75) -- 
#     ##    may need additional tuning; e.g. for eco (k_clust=450), k_clust/2 is too high
#     (clusters, n_improvements, n_tries, scores) = 
    floc_update(clusters, iter, all_genes, ratios, string_net) ##, max_improvements_per_iter)
#     (weight_r, weight_n, weight_m, weight_c, weight_v, weight_g) = get_score_weights( iter )
#     (weight_r_new, weight_n_new, weight_m_new, weight_c_new, weight_v_new, weight_g_new) = get_score_weights( iter + 1 )
#     n_motifs = get_n_motifs()
#     n_motifs_new = get_n_motifs(iter+1, n_iters)
#     if ( abs(weight_n) <= 0 && abs(weight_n_new) > 0 ) || ( weight_m <= 0 && weight_m_new > 0 ) || ( n_motifs != n_motifs_new )
#         for k in 1:k_clust clusters[k].changed[1] = true; end ## in this instance, need to force update of all clusters
#     end
#     iter += 1 ## now update clusters as if the new iteration has started
#     changed_rows = sum( [clusters[k].changed[1] for k=1:k_clust] )
#     changed_cols = sum( [clusters[k].changed[2] for k=1:k_clust] )
#     ## First, do the meme/mast-ing in parallel (only if m0 > 0)
#     (weight_r, weight_n, weight_m, weight_c, weight_v) = get_score_weights()
#     clusters = re_seed_all_clusters_if_necessary(clusters) ## avoid meme-ing 0-gene clusters
#     if weight_m > 0 
#         if nprocs() <= 1 clusters = re_meme_all_biclusters(clusters, false)
#         else clusters = re_meme_all_biclusters_parallel(clusters, false); end
#     end
#     ## Next fill the clusters' scores (in parallel)
#     if nprocs() <= 1 clusters = fill_all_cluster_scores( clusters, false, false );
#     else clusters = fill_all_cluster_scores_parallel( clusters, false, false ); end
#     println( "ITER: ", iter )
#     println( @sprintf( "r0: %.3f; n0: %.3f; m0: %.3f; c0: %.3f; v0: %.3f, g0: %.3f", weight_r, weight_n, weight_m, 
#                       weight_c, weight_v, weight_g ) )
#     println( "N_MOVES: ", n_tries )
#     println( "N_IMPROVEMENTS: ", n_improvements )
#     stats_df = print_cluster_stats(clusters)
#     stats_df["N_MOVES"] = n_tries
#     stats_df["N_IMPROVEMENTS"] = n_improvements
#     stats_df["N_CLUSTS_CHANGED_ROWS"] = changed_rows
#     stats_df["N_CLUSTS_CHANGED_COLS"] = changed_cols
#     println( "N_CLUSTS_CHANGED (ROWS): ", stats_df["N_CLUSTS_CHANGED_ROWS"] )
#     println( "N_CLUSTS_CHANGED (COLS): ", stats_df["N_CLUSTS_CHANGED_COLS"] )
#     gc()
#     (clusters, n_improvements, stats_df)
# end
