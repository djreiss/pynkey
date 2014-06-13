# Use FLOC algorithm to update clusters

import pandas as pd
import numpy as np
from numpy import random as rnd
from numpy import nan as NA

import utils as ut
import scores
from Bicluster import bicluster ##,fill_all_cluster_scores_par
import Bicluster as bic
import funcs
import meme
from params import n_iters,nthreads
import globals as glb

print 'importing floc'

## Get gain scores for all possible row/col moves
## Default: allow best 5 row- and 20 col- moves per bicluster instead of all of them!
def get_scores_all(clusters, iter, all_genes, ratios, string_net, max_row=9999, max_col=9999): ##5, max_col=20):
    ## First, collate move scores into a single DataFrame for stochastic sorting
    ## Use pd.DataFrame for scores rather than matrix

    ##bicluster.fill_all_cluster_scores(clusters, all_genes, ratios, string_net, ratios.columns.values)
    clusters = bic.fill_all_cluster_scores_par(clusters, threads=nthreads)

    counts_g = bicluster.get_all_cluster_row_counts( clusters, all_genes ) ## Clusters per gene counts - precompute
    tmp = [ clusters[k].get_floc_scoresDF_rows(iter, all_genes, ratios, counts_g) for k in clusters ]
    scoresDF_r = pd.concat(tmp, ignore_index=True)
    del tmp

    ## STANDARDIZE ROWS AND COLUMNS SEPARATELY!!!   -- is this really what we want?
    tmp1 = scores.get_combined_scores( scoresDF_r, iter, ratios )
    tmp1 = ut.sdize_vector( tmp1 )
    scoresDF_r['combined'] = pd.Series(tmp1)

    tmp = [ clusters[k].get_floc_scoresDF_cols(iter, ratios) for k in clusters ]
    scoresDF_c = pd.concat(tmp, ignore_index=True) ## This is much faster than mapreduce()
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
        scoresDF_r = pd.concat( shrunk_r, ignore_index=True )

    if max_col < 9999:
        tmp_c = scoresDF_c.groupby( 'k' )
        tmp_c = dict(list(tmp_c))
        shrunk_c = [pd.DataFrame() for i in range(len(tmp_c))]
        for i in range(len(tmp_c)):
            shrunk_c[i] = tmp_c[i].sort( 'combined' )[ 0:max_col ]
        scoresDF_c = pd.concat( shrunk_c, ignore_index=True )

    all_scores = pd.concat([scoresDF_r, scoresDF_c], ignore_index=True)
    return all_scores

# ## Instead of scores for ALL possible moves, make a matrix of n_best best scores for each row/col
def get_scores_best( all_scores, n_best_row=3, n_best_col=3 ):
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

    scores_out = pd.concat( [ pd.concat( dfs_r ), pd.concat( dfs_c ) ], ignore_index=True )
    return scores_out

## Native code: 119s
def rnd_bubblesort( scores, Nrepeats ): ##=None ): ## make sure scores is a copy, b/c NaNs will get replaced in this copy
    lsc = len(scores)
    if Nrepeats == None: ## is None:
        Nrepeats = lsc * 2
    ords = np.arange(lsc)
    rnd.shuffle(ords) ## start w/ random order
    (tmp1, tmp2) = (np.nanmin(scores), np.max(scores)) ##utils.minmax(scores)
    R = 2.0 * ( tmp2-tmp1 ) ## Denominator of value to compute prob. from
    the_max = tmp2
    n = lsc - 1
    sc = scores.copy()
    sc[ np.isnan(sc) ] = the_max ## replace NaN with maximum score
    n_switches = 0
    o1 = o2 = 0
    g1 = g2 = p = 0.
    for i in xrange(Nrepeats): ## TBD: this double loop can be sped up with weave!!!
        rnds = rnd.rand(n)
        for j in xrange(n):
            o1 = ords[j]
            o2 = ords[j+1]
            g1 = sc[o1]
            g2 = sc[o2]
            if g1 == g2 and g2 == the_max:
                continue
            p = 0.5 + ( g1 - g2 ) / R ## compute prob of switching
            if p > rnds[j]: ##rnd.rand() < p: ## switch???
                ords[j] = o2
                ords[j+1] = o1
                n_switches += 1
        if i % 1000 == 1:
            print( i, n_switches, Nrepeats )
    return ords

## This is jit-c'd via LVVM using numba - takes about 127 secs
## seems that it is not really jit-ed
## seems that I need to do more trickery to get it to compile; see
## http://www.reddit.com/r/Python/comments/1wnt5h/numba_code_slower_than_pure_python_code/
## https://stackoverflow.com/questions/18071285/numba-slow-when-assigning-to-an-array
## https://stackoverflow.com/questions/21468170/numba-code-slower-than-pure-python
#from numba import float64, int64, jit
#rnd_bubblesort2 = jit( int64[:]( float64[:], int64 ), nopython=True )( rnd_bubblesort )

## this is weaved (using c++) -- takes about 3.04 seconds
import scipy.weave
def rnd_bubblesort3( scores, Nrepeats=None ): ## make sure scores is a copy, b/c NaNs will get replaced in this copy
    lsc = len(scores)
    if Nrepeats is None:
        Nrepeats = lsc * 2
    ords = np.arange(lsc)
    rnd.shuffle(ords) ## start w/ random order
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
      int o1 = ords(j), o2 = ords(j+1);
      double g1 = sc(o1), g2 = sc(o2);
      if ( g1 == g2 && g2 == (double)the_max ) continue;
      double p = 0.5 + ( g1 - g2 ) / (double)R; // compute prob of switching
      if ( rnds(j) < p ) { // switch???
          ords(j) = o2;
          ords(j+1) = o1;
          switchesN(0) ++;
      }
    }
 """
        scipy.weave.inline(code,['ords','n','sc','R','switchesN','rnds','the_max'],
                           type_converters=scipy.weave.converters.blitz)
        if i % 1000 == 1:
            print i, switchesN[0], Nrepeats
    return ords

# ## TODO: add max_improvements param (to prevent really fast optimization at beginning before motifing turns on)
def update(clusters, iter, all_genes, ratios, string_net, max_no_improvements=250):

    scores_all = get_scores_all(clusters, iter, all_genes, ratios, string_net)

    scores_all2 = get_scores_best(scores_all)

    ## For FLOC, identify the BEST score for each row/column, then bubble sort those to get the order in
    ##    which to perform the moves.
    ## Note this is wrong right now - it sorts ALL k scores for each row/col. 
    ##  Need to just use the BEST score for each row/col and then bubblesort these.
    ord = rnd_bubblesort3( scores_all2['combined'].values ) ##, n_sort_iter)
    ##print scores_all2.ix[ord,:].head(); print scores_all2.ix[ord,:].tail()

    new_clusters = saved_clusters = funcs.copy_clusters( clusters, deep=True ) ## make a copy for updating
    ## don't need to update n or m scores if their weights are 0:
    weight_r, weight_n, weight_m, weight_c, weight_v, weight_g = scores.get_score_weights( iter, ratios )

    all_resids = np.array( [ clusters[k].resid for k in clusters.keys() ] )
    all_dens_string = np.array( [ clusters[k].dens_string for k in clusters.keys() ] )
    all_meanp_meme = np.array( [ clusters[k].meanp_meme for k in clusters.keys() ] )
    all_volumes = np.zeros_like( all_resids )
    mean_resid = np.nanmean(all_resids)
    mean_dens_string = np.nanmean(all_dens_string) ## Should usually increase, sometimes get worse
    mean_meanp_meme = np.nanmean(all_meanp_meme) ## Should usually decrease, sometimes get worse
    counts_gene = bicluster.get_all_cluster_row_counts( clusters, all_genes ) ## Clusters per gene score
    #scores_gene = get_cluster_row_count_scores( counts_gene )
    #mean_score_gene = nanmean(scores_gene)
    no_improvements = 0 ## Number of updates in which we haven't seen an improvement
    n_improvements = n_tries = 0

    for i in ord: ## Update bicusters -- should "store" the one with the best mean resid during the updates
        sc = scores_all2.ix[i]
        kUpd = sc.k
        cc = new_clusters[kUpd].copy(deep=False)
        row_col = sc.row_col
        if sc.is_row_col == 'r':
            ##if sc.is_in && len(cc.rows) <= min_rows: continue ## Don't remove if we're already at the min. Now: use volume score instead
            cc.rows = cc.rows[ cc.rows != row_col ] if sc.is_in else np.append( cc.rows, row_col )
            cc.changed[0] = True
        elif sc.is_row_col == 'c':
            cc.cols = cc.cols[ cc.cols != row_col ] if sc.is_in else np.append( cc.cols, row_col )
            cc.changed[1] = True

        cc.resid = cc.compute_residue( ratios )
        all_resids[kUpd] = cc.resid
        if sc.is_row_col == 'r':  ## Only update network/motif scores if it's a row, duh!
            if abs(weight_n) > 0:            ## only compute if need to
                cc.dens_string = cc.compute_network_density( string_net )
                all_dens_string[kUpd] = cc.dens_string
            if weight_m > 0:            ## only compute if need to
                cc.meanp_meme = cc.compute_meme_pval()
                all_meanp_meme[kUpd] = cc.meanp_meme
            ## TODO: incorporate score_g (nclust per gene) and score_v (ngene per clust) into combined scores
#            if weight_g > 0:
#                counts_gene[ row_col ] += (-1 if sc.is_in else 1) ## Clusters per gene score
#                mean_score_gene = nanmean(get_cluster_row_count_scores(cc, counts_gene)) ## Clusters per gene score
#                counts_gene[ row_col ] -= (-1 if scis_in else 1) ## revert it back
        ##all_scores[kUpd] = get_combined_score( cc.resid, cc.dens_string, cc.meanp_meme, 0.0f0 )
        ##all_scores = get_combined_scores( sdize_vector( all_resids ), sdize_vector( all_dens_string ), 
        ##                                 sdize_vector( all_meanp_meme ), all_volumes )
        ## Need to normalize these deltas vs. the stddev of the corresponding scores in the "score" dataframe!
        score_delta = scores.get_combined_score( (np.nanmean(all_resids) - mean_resid) / np.nanstd(all_resids), 
                                   (np.nanmean(all_dens_string) - mean_dens_string) / np.nanstd(all_dens_string),
                                   (np.nanmean(all_meanp_meme) - mean_meanp_meme) / np.nanstd(all_meanp_meme), 
                                   0., 0., ##(np.nanmean(scores_gene) - mean_score_gene) / np.nansd(scores_gene)
                                                 iter, ratios )
        new_clusters[kUpd] = cc
        n_tries += 1
        ##print np.nanmean(all_scores), mean_scores
        if score_delta <= 0: ## mean bicluster scores improved, so store the entire cluster stack
            n_improvements += 1
            mean_resid = np.nanmean(all_resids) ## Should usually decrease, sometimes get worse
            mean_dens_string = np.nanmean(all_dens_string) ## Should usually increase, sometimes get worse
            mean_meanp_meme = np.nanmean(all_meanp_meme) ## Should usually decrease, sometimes get worse
            ##mean_scores = np.nanmean(all_scores) ## Should usually decrease, sometimes get worse
            ##saved_clusters[kUpd] = cc ##copy_cluster(cc, true)
            saved_clusters = funcs.copy_clusters( new_clusters, deep=False ) ## make a copy to keep the best update
            output = '%d %.4f %.4f %.4f %.4f %d %c %s %s %.4f %.4f %.4f %d %d %d' % \
                (i, mean_resid, mean_dens_string, mean_meanp_meme, score_delta, ##mean_scores,
                 kUpd, sc.is_row_col, 'remove' if sc.is_in else 'add', row_col, 
                 cc.resid, clusters[kUpd].resid, cc.resid-clusters[kUpd].resid, 
                 n_tries, n_improvements, no_improvements)
            print output
            no_improvements = 0
        else:
            no_improvements += 1
        if no_improvements > max_no_improvements: 
            break
    return (saved_clusters, n_improvements, n_tries, scores_all2.ix[ord,:])

def run(clusters, iter, all_genes, ratios, string_net):
    ##clusters = fill_cluster_scores(clusters) ## dont need this since each clust's scores are updated in floc.update
    ## allow more updates if there are more clusters??? Tuned to k_clust/2 for Hpy (where k_clust is 75) -- 
    ##    may need additional tuning; e.g. for eco (k_clust=450), k_clust/2 is too high
    clusters, n_improvements, n_tries, scores_all = update(clusters, iter, all_genes, ratios, string_net)
    glb.clusters = clusters

    weight_r, weight_n, weight_m, weight_c, weight_v, weight_g = scores.get_score_weights( iter, ratios )
    weight_r_new, weight_n_new, weight_m_new, weight_c_new, weight_v_new, weight_g_new = \
        scores.get_score_weights( iter + 1, ratios )
    n_motifs = meme.get_n_motifs(iter, n_iters)
    n_motifs_new = meme.get_n_motifs(iter+1, n_iters)
    if ( abs(weight_n) <= 0 and abs(weight_n_new) > 0 ) or \
            ( weight_m <= 0 and weight_m_new > 0 ) or \
            ( n_motifs != n_motifs_new ):
        for k in clusters.keys():
            clusters[k].changed[0] = True ## in this instance, need to force update of all clusters
    iter += 1 ## now update clusters as if the new iteration has started
    changed_rows = sum( [clusters[k].changed[0] for k in range(len(clusters))] )
    changed_cols = sum( [clusters[k].changed[1] for k in range(len(clusters))] )

    ## First, do the meme/mast-ing in parallel (only if m0 > 0)
    ## avoid meme-ing 0-gene clusters
    clusters = funcs.re_seed_all_clusters_if_necessary(clusters, ratios, all_genes, min_rows=3, max_rows=80 )
    glb.clusters = clusters ## need to do this before the parallelized call:
    clusters = bic.fill_all_cluster_scores_par(clusters, threads=nthreads)
    glb.clusters = clusters
    print 'ITER:', iter
    print 'r0: %.3f; n0: %.3f; m0: %.3f; c0: %.3f; v0: %.3f, g0: %.3f' % ( weight_r, weight_n, weight_m, 
                                                                           weight_c, weight_v, weight_g )
    print 'N_MOVES:', n_tries
    print 'N_IMPROVEMENTS:', n_improvements
    stats_df = funcs.print_cluster_stats(clusters, ratios, iter, glb.startTime )
    stats_df['N_MOVES'] = n_tries
    stats_df['N_IMPROVEMENTS'] = n_improvements
    stats_df['N_CLUSTS_CHANGED_ROWS'] = changed_rows
    stats_df['N_CLUSTS_CHANGED_COLS'] = changed_cols
    print 'N_CLUSTS_CHANGED (ROWS): ', stats_df.N_CLUSTS_CHANGED_ROWS[0]
    print 'N_CLUSTS_CHANGED (COLS): ', stats_df.N_CLUSTS_CHANGED_COLS[0]
    return (clusters, n_improvements, stats_df)
