# Use FLOC algorithm to update clusters

function get_floc_scoresDF_rows(cluster::bicluster, counts_g::Vector{Int32})
    inds = [1:length(cluster.scores_r)]
    is_in_r = in( inds, cluster.rows )
    (weight_r, weight_n, weight_m, weight_c, weight_v, weight_g) = get_score_weights() ## only use them if their weights are > 0
    NAs = weight_r <= 0 || weight_n <= 0 || weight_m <= 0 ? fill(NA, length(cluster.scores_r)) : NA
    score_r = weight_r > 0 ? cluster.scores_r : NAs
    score_n = abs(weight_n) > 0 ? cluster.scores_n : NAs
    score_m = weight_m > 0 ? cluster.scores_m : NAs
    score_vr = get_cluster_volume_row_scores( cluster, is_in_r )
    score_g = get_cluster_row_count_scores( cluster, counts_g, is_in_r )
    out = DataFrame( { "row_col_ind" => inds,
                      "is_in" => is_in_r,
                      "is_row_col" => fill('r', length(cluster.scores_r)), ## CANT: move this outside the loop
                      "k" => fill(cluster.k, length(cluster.scores_r)),
                      "score" => score_r,
                      "score_n" => score_n,
                      "score_m" => score_m,
                      "score_v" => score_vr,
                      "score_g" => score_g } )
    out
end

function get_floc_scoresDF_cols(cluster::bicluster)
    inds = [1:length(cluster.scores_c)]
    is_in_c = in( inds, cluster.cols )
    (weight_r, weight_n, weight_m, weight_c, weight_v, weight_g) = get_score_weights() ## only use them if their weights are > 0
    NAs = fill(NA, length(cluster.scores_c))
    score_c = weight_c > 0 ? cluster.scores_c : NAs
    score_vc = get_cluster_volume_col_scores( cluster, is_in_c )
    out = DataFrame( { "row_col_ind" => inds,
                      "is_in" => is_in_c,
                      "is_row_col" => fill('c', length(cluster.scores_c)), ## CANT: move this outside the loop
                      "k" => fill(cluster.k, length(cluster.scores_c)),
                      "score" => score_c,
                      "score_n" => NAs, ## CANT: move this outside the loop
                      "score_m" => NAs,  ## CANT: move this outside the loop
                      "score_v" => score_vc,
                      "score_g" => NAs } )
    out
end

## Default: allow best 5 row- and 20 col- moves per bicluster instead of all of them!
#get_floc_scores_all(clusters::Dict{Int64,bicluster}) = get_floc_scores_all(clusters, 9999, 9999) ## 5, 20)

## Get gain scores for all possible row/col moves
function get_floc_scores_all(clusters::Dict{Int64,bicluster}, max_row::Int64=9999, max_col::Int64=9999) 
    global k_clust

    ## First, collate move scores into a single DataFrame for stochastic sorting
    ## Use DataFrames for scores rather than matrix

    counts_g::Vector{Int32} = get_cluster_row_counts( clusters ) ## Clusters per gene counts - precompute

    #scoresDF_r::DataFrame = mapreduce( k->get_floc_scoresDF_rows(clusters[k], k), vcat, 1:k_clust )
    tmp::Vector{DataFrame} = [get_floc_scoresDF_rows(clusters[k], counts_g) for k=1:k_clust]
    scoresDF_r::DataFrame = rbind(tmp) ## This is much faster than mapreduce() above!!!

    ## NEED TO STANDARDIZE ROWS AND COLUMNS SEPARATELY!!!   -- is this really what we want?
    tmp_r = DataFrame( { "score" => convert(Vector{Float32}, scoresDF_r["score"] ),
                        "score_n" => convert(Vector{Float32}, scoresDF_r["score_n"] ), 
                        "score_m" => convert(Vector{Float32}, scoresDF_r["score_m"] ),
                        "score_v" => convert(Vector{Float32}, scoresDF_r["score_v"] ),
                        "score_g" => convert(Vector{Float32}, scoresDF_r["score_g"] ) } )
    ## Weights are set here: 
    tmp1 = get_combined_scores( sdize_vector( tmp_r["score"].data ), sdize_vector( tmp_r["score_n"].data ), 
                               sdize_vector( tmp_r["score_m"].data ), ##sdize_vector( tmp_r["score_v"].data ),
                               ##sdize_vector( scores_g ) )
                               sdize_vector(tmp_r["score_v"].data),
                               sdize_vector(tmp_r["score_g"].data))
    tmp1 = sdize_vector( tmp1 )
    scoresDF_r["combined"] = tmp1

    ##scoresDF_c::DataFrame = mapreduce( k->get_floc_scoresDF_cols(clusters[k], k), vcat, 1:k_clust )
    tmp = [get_floc_scoresDF_cols(clusters[k]) for k=1:k_clust]
    scoresDF_c::DataFrame = rbind(tmp) ## This is much faster than mapreduce() above!!!
    tmp_na::Vector{Float32} = float32( fill(NA, size(scoresDF_c,1)) )

    tmp_c = DataFrame( { "score" => convert(Vector{Float32}, scoresDF_c["score"] ),
                        "score_n" => tmp_na,
                        "score_m" => tmp_na,
                        "score_v" => convert(Vector{Float32}, scoresDF_c["score_v"] ),
                        "score_g" => tmp_na } )
    ## Weights are set here:
    tmp1 = get_combined_scores( sdize_vector( tmp_c["score"].data ), tmp_na, ##tmp_c["score_n"].data, 
                               tmp_na, ##tmp_c["score_m"].data, ##sdize_vector( tmp_c["score_v"].data ), 
                               sdize_vector(tmp_c["score_v"].data), tmp_na ) ##tmp_c["score_g"].data )
    (weight_r, weight_n, weight_m, weight_c, weight_v, weight_g) = get_score_weights() ## get weight_c to up-weight cols. vs. rows, if desired
    tmp1 = sdize_vector( tmp1 ) .* weight_c
    scoresDF_c["combined"] = tmp1

    if max_row <= maximum(scoresDF_r["row_col_ind"])
        tmp_r = groupby( scoresDF_r, "k" )
        shrunk_r = Array(DataFrame, length(tmp_r))
        for i in 1:length(tmp_r) 
            ord = sortperm( tmp_r[i]["combined"] )[1:max_row]
            shrunk_r[i] = tmp_r[i][ ord, : ]
        end
        scoresDF_r = rbind( shrunk_r )
    end

    if max_row <= maximum(scoresDF_r["row_col_ind"])
        tmp_c = groupby( scoresDF_c, "k" )
        shrunk_c = Array(DataFrame, length(tmp_c))
        for i in 1:length(tmp_c) 
            ord = sortperm( tmp_c[i]["combined"] )[1:max_col]
            shrunk_c[i] = tmp_c[i][ ord, : ]
        end
        scoresDF_c = rbind( shrunk_c )
    end

    scores = vcat(scoresDF_r, scoresDF_c)
    #println(scores)
    #write_table("output/$(organism)_scores.tsv",scores) 
    scores
end

#get_floc_scores_best( scores::DataFrame ) = get_floc_scores_best( scores, 3, 3 ) ## default n_best=2

## Instead of scores for ALL possible moves, make a matrix of n_best best scores for each row/col
function get_floc_scores_best( scores::DataFrame, n_best_row::Int64=3, n_best_col::Int64=3 )
    dfs_r::Vector{AbstractDataFrame} = [ sub( scores, :(is_row_col .== 'r') ) ]
    if n_best_row < 9999  ## Now uses the nice functions in DataFrame like with() or group() or some such
        tmp = groupby( dfs_r[1], "row_col_ind" ) ## Cool!
        dfs_r = Array(AbstractDataFrame, length(tmp))
        for r in 1:length(tmp)
            tmp2 = tmp[r] ##sub( tmp[r], :(row_col_ind .== row_col_ind[1]) )
            if n_best_row == 1 ord = findmin( tmp2["combined"] )[2]
            elseif n_best_row >= size(tmp2, 1) ord = 1:size(tmp2, 1)
            else ord = sortperm( tmp2["combined"] )[ 1:n_best_row ]
            end
            dfs_r[r] = tmp2[ ord, : ]
        end
    end

    dfs_c::Vector{AbstractDataFrame} = [ sub( scores, :(is_row_col .== 'c') ) ]
    if n_best_col < 9999  ## Now uses the nice functions in DataFrame like with() or group() or some such
        tmp = groupby( dfs_c[1], "row_col_ind" ) ## Cool!
        dfs_c = Array(AbstractDataFrame, length(tmp))
        for c in 1:length(tmp)
            tmp2 = tmp[c] ##sub( tmp[c], :(row_col_ind .== $c) )
            if n_best_col == 1 ord = findmin( tmp2["combined"] )[2]
            elseif n_best_col >= size(tmp2, 1) ord = 1:size(tmp2, 1)
            else ord = sortperm( tmp2["combined"] )[ 1:n_best_col ]
            end
            dfs_c[c] = tmp2[ ord,: ]
        end
    end
    scores2 = rbind( rbind( dfs_r ), rbind( dfs_c ) )
    scores2
end

#rnd_bubblesort( scores::Vector{Float32} ) = rnd_bubblesort( scores, nothing )

function rnd_bubblesort( scores::Vector{Float32}, Nrepeats=nothing )
    const lsc = length(scores)
    if Nrepeats == nothing Nrepeats = lsc * 2; end
    ord::Vector{Int32} = shuffle!( [1:lsc] ) ## start w/ random order
    tmp = range(scores); const R=2.0 * ( tmp[2]-tmp[1] ) ## Denominator of value to compute prob. from
    the_max::Float32 = tmp[2]
    const n = lsc - 1
    n_switches::Int = 0
    o1::Int32 = 0; o2::Int32 = 0; g1::Float32 = 0.; g2::Float32 = 0.; p::Float32 = 0.
    for i=1:Nrepeats
        n_switches = 0
        for j=1:n
            o1 = ord[j]
            o2 = ord[j+1]
            g1 = scores[o1]; if g1 != g1 g1 = the_max; end ## replace NaN with maximum score
            g2 = scores[o2]; if g2 != g2 g2 = the_max; end ## replace NaN with maximum score
            if g1 == g2 && g2 == the_max continue; end
            p = 0.5 + ( g1 - g2 ) / R ## compute prob of switching
            if rand() < p ## switch???
                ord[j] = o2
                ord[j+1] = o1
                n_switches += 1
            end
        end
        if i % 10000 == 1 println(i," ",n_switches," ",lsc*2); end
    end
    ord
end

#floc_update(clusters) = floc_update(clusters, 25)

## TODO: add max_improvements param (to prevent really fast optimization at beginning before motifing turns on)
function floc_update(clusters::Dict{Int64,bicluster}, max_no_improvements=25)
    global k_clust, ratios, string_net

    scores::DataFrame = get_floc_scores_all(clusters);
    scores2::DataFrame = get_floc_scores_best(scores);

    ## For FLOC, identify the BEST score for each row/column, then bubble sort those to get the order in
    ##    which to perform the moves.
    ## Note this is wrong right now - it sorts ALL k scores for each row/col. 
    ##  Need to just use the BEST score for each row/col and then bubblesort these.
    ord::Vector{Int32} = rnd_bubblesort(convert(Vector{Float32}, scores2["combined"])) ##, n_sort_iter) 
    println( head(scores2[ord,:]), "\n", tail(scores2[ord,:]) )

    new_clusters = saved_clusters = copy_clusters( clusters, true, false ); ## make a copy for updating
    (weight_r, weight_n, weight_m, weight_c, weight_v, weight_g) = get_score_weights() ## DONE: don't need to update n or m scores if their weights are 0

    all_resids = float32( [ clusters[k].resid for k=1:k_clust ] )
    all_dens_string = float32( [ clusters[k].dens_string for k=1:k_clust ] )
    all_meanp_meme = float32( [ clusters[k].meanp_meme for k=1:k_clust ] )
    all_volumes = zeros( Float32, length(all_resids) )
    mean_resid = nanmean(all_resids)
    mean_dens_string = nanmean(all_dens_string) ## Should usually increase, sometimes get worse
    mean_meanp_meme = nanmean(all_meanp_meme) ## Should usually decrease, sometimes get worse
    counts_gene = get_cluster_row_counts( clusters ) ## Clusters per gene score
    #scores_gene = get_cluster_row_count_scores( counts_gene )
    #mean_score_gene = nanmean(scores_gene)
    no_improvements = 0 ## Number of updates in which we haven't seen an improvement
    n_improvements = 0; n_tries = 0
    for i=ord ## Update bicusters -- should "store" the one with the best mean resid during the updates
        sc = scores2[i,:]
        kUpd = sc["k"][1]
        cc::bicluster = copy_cluster(new_clusters[kUpd], true, true)
        row_col = sc["row_col_ind"][1]
        if sc["is_row_col"][1] == 'r'
            ##if sc["is_in"][1] && length(cc.rows) <= min_rows continue; end ## Don't remove if we're already at the min. Now: use volume score instead
            cc.rows = sc["is_in"][1] ? remove(cc.rows, row_col) : [cc.rows, row_col]
            cc.changed[1] = true
        elseif sc["is_row_col"][1] == 'c'
            cc.cols = sc["is_in"][1] ? remove(cc.cols, row_col) : [cc.cols, row_col]
            cc.changed[2] = true
        end
        cc.resid = bicluster_residue( cc, ratios )
        all_resids[kUpd] = cc.resid
        if sc["is_row_col"][1] == 'r'  ## Only update network/motif scores if it's a row, duh!
            if abs(weight_n) > 0            ## only compute if need to
                cc.dens_string = bicluster_network_density( cc, string_net )
                all_dens_string[kUpd] = cc.dens_string
            end
            if weight_m > 0            ## only compute if need to
                cc.meanp_meme = bicluster_meme_pval( cc )
                all_meanp_meme[kUpd] = cc.meanp_meme
            end
            ## TODO: incorporate score_g (nclust per gene) and score_v (ngene per clust) into combined scores
#             if weight_g > 0 
#                 counts_gene[ row_col ] += (sc["is_in"][1] ? -1 : 1) ## Clusters per gene score
#                 mean_score_gene = nanmean(get_cluster_row_count_scores(cc, counts_gene)) ## Clusters per gene score
#                 counts_gene[ row_col ] -= (sc["is_in"][1] ? -1 : 1) ## revert it back
#             end
        end
        ##all_scores[kUpd] = get_combined_score( cc.resid, cc.dens_string, cc.meanp_meme, 0.0f0 )
        ##all_scores = get_combined_scores( sdize_vector( all_resids ), sdize_vector( all_dens_string ), 
        ##                                 sdize_vector( all_meanp_meme ), all_volumes )
        ## Need to normalize these deltas vs. the stddev of the corresponding scores in the "score" dataframe!
        score_delta = get_combined_score( (nanmean(all_resids) - mean_resid) / nansd(all_resids), 
                                         (nanmean(all_dens_string) - mean_dens_string) / nansd(all_dens_string),
                                         (nanmean(all_meanp_meme) - mean_meanp_meme) / nansd(all_meanp_meme), 
                                         0.0f0, 0.0f0 ) ##(nanmean(scores_gene) - mean_score_gene) / nansd(scores_gene) )
        new_clusters[kUpd] = cc
        n_tries += 1
        ##println(nanmean(all_scores)," ",mean_scores)
        if score_delta <= 0 ## mean bicluster scores improved, so store the entire cluster stack
            n_improvements += 1
            mean_resid = nanmean(all_resids) ## Should usually decrease, sometimes get worse
            mean_dens_string = nanmean(all_dens_string) ## Should usually increase, sometimes get worse
            mean_meanp_meme = nanmean(all_meanp_meme) ## Should usually decrease, sometimes get worse
            ##mean_scores = nanmean(all_scores) ## Should usually decrease, sometimes get worse
            ##saved_clusters[kUpd] = cc ##copy_cluster(cc, true)
            saved_clusters = copy_clusters( new_clusters, false, false ) ## make a copy to keep the best update
            output = @sprintf( "%d %.4f %.4f %.4f %.4f %d %c %s %d %.4f %.4f %.4f %d %d %d",
                              i, mean_resid, mean_dens_string, mean_meanp_meme, score_delta, ##mean_scores,
                              Base.int(kUpd), sc["is_row_col"][1], sc["is_in"][1]?"remove":"add", row_col, 
                              cc.resid, clusters[kUpd].resid, cc.resid-clusters[kUpd].resid, 
                              n_tries, n_improvements, no_improvements)
            println( output )
            no_improvements = 0
        else
            no_improvements += 1
        end
        if no_improvements > max_no_improvements break; end
    end
    (saved_clusters, n_improvements, n_tries, scores2[ord,:])
end

function do_floc(clusters::Dict{Int64,bicluster})
    global iter, k_clust, max_improvements_per_iter
    ##clusters = fill_cluster_scores(clusters) ## dont need this since each cluster's scores are updated in floc_update
    ## allow more updates if there are more clusters??? Tuned to k_clust/2 for Hpy (where k_clust is 75) -- 
    ##    may need additional tuning; e.g. for eco (k_clust=450), k_clust/2 is too high
    (clusters, n_improvements, n_tries, scores) = floc_update(clusters, max_improvements_per_iter); 
    (weight_r, weight_n, weight_m, weight_c, weight_v, weight_g) = get_score_weights( iter )
    (weight_r_new, weight_n_new, weight_m_new, weight_c_new, weight_v_new, weight_g_new) = get_score_weights( iter + 1 )
    n_motifs = get_n_motifs()
    n_motifs_new = get_n_motifs(iter+1, n_iters)
    if ( abs(weight_n) <= 0 && abs(weight_n_new) > 0 ) || ( weight_m <= 0 && weight_m_new > 0 ) || ( n_motifs != n_motifs_new )
        for k in 1:k_clust clusters[k].changed[1] = true; end ## in this instance, need to force update of all clusters
    end
    iter += 1 ## now update clusters as if the new iteration has started
    changed_rows = sum( [clusters[k].changed[1] for k=1:k_clust] )
    changed_cols = sum( [clusters[k].changed[2] for k=1:k_clust] )
    ## First, do the meme/mast-ing in parallel (only if m0 > 0)
    (weight_r, weight_n, weight_m, weight_c, weight_v) = get_score_weights()
    clusters = re_seed_all_clusters_if_necessary(clusters) ## avoid meme-ing 0-gene clusters
    if weight_m > 0 
        if nprocs() <= 1 clusters = re_meme_all_biclusters(clusters, false)
        else clusters = re_meme_all_biclusters_parallel(clusters, false); end
    end
    ## Next fill the clusters' scores (in parallel)
    if nprocs() <= 1 clusters = fill_all_cluster_scores( clusters, false, false );
    else clusters = fill_all_cluster_scores_parallel( clusters, false, false ); end
    println( "ITER: ", iter )
    println( @sprintf( "r0: %.3f; n0: %.3f; m0: %.3f; c0: %.3f; v0: %.3f, g0: %.3f", weight_r, weight_n, weight_m, 
                      weight_c, weight_v, weight_g ) )
    println( "N_MOVES: ", n_tries )
    println( "N_IMPROVEMENTS: ", n_improvements )
    stats_df = print_cluster_stats(clusters)
    stats_df["N_MOVES"] = n_tries
    stats_df["N_IMPROVEMENTS"] = n_improvements
    stats_df["N_CLUSTS_CHANGED_ROWS"] = changed_rows
    stats_df["N_CLUSTS_CHANGED_COLS"] = changed_cols
    println( "N_CLUSTS_CHANGED (ROWS): ", stats_df["N_CLUSTS_CHANGED_ROWS"] )
    println( "N_CLUSTS_CHANGED (COLS): ", stats_df["N_CLUSTS_CHANGED_COLS"] )
    gc()
    (clusters, n_improvements, stats_df)
end
