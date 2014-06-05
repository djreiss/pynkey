import warnings

import numpy as np

def copy_clusters( clusters, deep=True ):
    new_clusters = {} ## make a copy for updating
    for k in clusters.keys():
        new_clusters[k] = clusters[k].copy( deep )
    return new_clusters

def re_seed_all_clusters_if_necessary( clusters, ratios ):
    for k in clusters.keys():
        clusters[k] = clusters[k].re_seed_if_necessary( clusters, ratios )
    return clusters

def matrix_residue( rats ):
    if np.ndim( rats ) < 2 or np.size( rats, 0 ) <= 1 or np.size( rats, 1 ) <= 1: ## or \
            ##np.mean( rats.isnull().values ) > 0.95:
        warnings.warn( "COULD NOT COMPUTE RESIDUE" )
        return 1.0

    rats = rats.values ## do it all in numpy - faster
    d_rows = np.nanmean(rats, 1) ##rats.mean(1)
    d_cols = np.nanmean(rats, 0) ##rats.mean(0)
    d_all = np.nanmean(d_rows) ##d_rows.mean() ## pandas default is to ignore nan's
    rats = rats + d_all - np.add.outer( d_rows, d_cols )
    ## another way of doing this, but about the same speed:
    ##rats = (rats+d_all-d_cols).T-d_rows

    average_r = np.nanmean( np.abs(rats) )
    return average_r

def matrix_var( rats, var_add=0.1 ):
    rats = rats.values
    mn = np.nanmean(rats, 0) ## subtract bicluster mean profile
    return np.nanvar(rats-mn) / (np.nanvar(mn) + var_add)

## Use edge density score: sum of edge weights (=number of edges for all weights =1) / number of nodes^2
## This is the faster version that uses SubDataFrames
## TODO: use numexpr to speed up and avoid temporary array creation?
def subnetwork_density( rows, network, already_subnetted=False ):
    ##net2 = net1[ np.in1d( net1.protein2, rows ) ]
    if already_subnetted: ## presumed already subnetted and index set to protein2
        net1 = network
    else:
        net1 = network.ix[rows]
        net1.set_index( ['protein2'], inplace=True ) ##[ np.in1d( net1.protein2, rows ) ]
    net2 = net1.ix[ rows ]
    dens = float(np.sum( net2.weight )) / (float(len(rows))**2) ## Already symmetrized, need to decrease count by 1/2
    return np.log10( dens+1e-9 )

# function fill_all_cluster_scores( clust::bicluster, iter::Int64, force=false, verbose=false )
#     global ratios, string_net
#     ##clust = re_seed_bicluster_if_necessary(clust)
#     if verbose println(clust.k, " ", length(clust.rows), " ", length(clust.cols)); end
#     (weight_r, weight_n, weight_m, weight_c, weight_v) = get_score_weights(iter)
#     #println("HERE5 $weight_r $(clust.changed)")
#     if force || ( ( weight_r + weight_c > 0 ) && ( sum(clust.changed) > 0 ) )  
#         try
#             clust = get_cluster_expr_rowcol_scores( clust, ratios ) ## Actually need to do it for rows too, even if only cols changed
#         catch x
#             warn( "ERROR WITH ROWCOL SCORES!" )
#             println( x )
#         end
#     end
#     #println("HERE6")
#     if force || clust.changed[1] ## rows only
#         if abs(weight_n) > 0 
#             try
#                 clust = get_cluster_network_row_scores( clust, string_net )
#             catch x
#                 warn( "ERROR WITH NETWORK SCORES!" )
#                 println( x )
#             end
#         end
#         #println("HERE6a")
#         if weight_m > 0 
#             try
#                 clust = get_cluster_meme_row_scores( clust )
#             catch x
#                 warn( "ERROR WITH MEME SCORES!" )
#                 println( x )
#             end
#         end
#     end 
#     #println("HERE7")
#     clust.changed[1] = clust.changed[2] = false
#     #println("HERE8")
#     clust
# end

# ## DONE: if anything needs to be parallelized, it's this function. (At least the meme-ing.) -- see re_meme_all_biclusters_parallel
# ## DONE: Also use the bicluster.changed field consistently and don't re-calc for biclusters that haven't changed.
# function fill_all_cluster_scores( clusters::Dict{Int64,bicluster}, force=false, verbose=false )
#     global k_clust
#     clusters = re_seed_all_clusters_if_necessary(clusters)    
#     ## Get all the scores for adding/removing rows/cols from each bicluster
#     for k=1:k_clust ## DONE: good possibility for parallelization? Note, we have to send over all the data
#         clust = clusters[k]     ## (including ratios, networks) to children.
#         if ! force && ! any(clust.changed) continue; end
#         if verbose println(k, " ", length(clust.rows), " ", length(clust.cols)); end
#         clust = fill_cluster_scores(clust, force)
#         clusters[ k ] = clust
#     end
#     clusters
# end

# function pre_load_child_nodes() 
#     global organism, k_clust, ratios, string_net, all_genes, n_iters
#     ## Need to pre-send all data to children -- I dont know how to send a variable in the @everywhere call so use a fixed filename
#     ## This only needs to be done once, doing it multiple times causes a big memory leak on the children.
#     println( "Sending data to child nodes" )
#     fname = "tmp_send_data.jldz"
#     save_jld( fname, (organism, k_clust, ratios, string_net, all_genes, n_iters, motif_width_range, distance_search, distance_scan) ) 
#                       #genome_seq, anno, op_table, allSeqs_fname, all_bgFreqs, all_genes, all_rows, iter, clusters) )
#     @everywhere (organism, k_clust, ratios, string_net, all_genes, n_iters, motif_width_range, distance_search, distance_scan) = load_jld("tmp_send_data.jldz")
#     sleep( 5 ) ## make sure children have time to load the data before the file is removed ???
#     rm( fname )
# end

# function fill_all_cluster_scores_parallel( clusters::Dict{Int64,bicluster}, force=false, verbose=false )
#     global k_clust, iter
#     clusters = re_seed_all_clusters_if_necessary(clusters) ## do this first since it relies on a global "clusters"
#     data::Array{Any,1} = []
#     for k in 1:k_clust
#         if ! force && ! any(clusters[k].changed) continue; end ## If bicluster not changed, skip it
#         b = bicluster( k, clusters[k].rows, clusters[k].cols ) ## copy the cluster but only send over the necessary parts
#         b.mast_out = clusters[k].mast_out
#         if ! force b.changed = clusters[k].changed
#         else b.changed[1] = b.changed[2] = true; end
#         dat = { "iter" => iter, "k" => b.k, "biclust" => b, "verbose" => verbose, "force" => force }
#         push!( data, dat )
#     end
#     new_clusts = pmap( fill_cluster_scores, data ) ## returns an Array{Any,1}
#     @everywhere gc() ## force a gc -- does it help?
#     for i in 1:length(new_clusts)
#         b = new_clusts[i]  ## new (updated) cluster
#         k = b.k
#         bb = clusters[k] ## original cluster
#         b.meme_out = bb.meme_out ## this did not get updated in b so just copy it over
#         if ! force ## THIS STUFF BELOW IS BECAUSE WE CALL fill_cluster_scores(x::Dict) with force=false
#             if ! any(bb.changed) ## rows AND cols not changed so scores_r, scores_c were NOT updated in b
#                 b.scores_r = bb.scores_r
#                 b.scores_c = bb.scores_c
#                 b.resid = bb.resid; b.var = bb.var
#             end
#             if ! bb.changed[1] ## rows not changed so scores_r, scores_n were NOT updated in b
#                 b.scores_m = bb.scores_m; b.meanp_meme = bb.meanp_meme
#                 b.scores_n = bb.scores_n; b.dens_string = bb.dens_string
#             end
#         end
#         clusters[k] = b 
#     end
#     clusters
# end

# function fill_cluster_scores( x::Dict{Any,Any} ) 
#     k = x["k"]
#     iter = x["iter"]
#     b = x["biclust"]
#     force = x["force"]
#     verbose = x["verbose"]
#     out = fill_cluster_scores(b, iter, force, verbose)   ### note changed force to false; this may be buggy
#     out
# end

## Find the flagellar cluster, whew!!!
# function clusters_w_func( func::ASCIIString, clusters, n_best=1 )
#     global anno, ratios
#     inds = findn([ismatch(Regex(func),anno["desc"][i]) for i=1:size(anno,1)] .== true)
#     r_rownames = rownames(ratios)
#     inds2 = int32([in(r_rownames, anno["sysName"][i]) ? 
#                    findn(r_rownames .== anno["sysName"][i])[1] : 0 for i=inds])
#     inds2 = inds2[ inds2 .!= 0 ]

#     ord = sortperm( int64([length(findin(clusters[k].rows,int64(inds2))) for k=1:length(clusters)]) )
#     kInds = ord[ (end-n_best+1):end ] ## The n_best clusters w/ the most genes annotated with "flagell"
#     ##kInd = findmax(int64([length(findin(clusters[k].rows,int64(inds2))) for k=1:length(clusters)]))[ 2 ]  
    
#     for kInd in kInds
#         genes = r_rownames[clusters[kInd].rows] ##rows]
#         ##println(genes) ## print the genes
#         genes = genes[ findin(genes, anno["sysName"].data) ]
#         println(kInd, "\n", anno[in(anno["sysName"].data,genes),["sysName","desc"]])
#     end
#     kInds
# end

# function print_cluster_stats( clusters::Dict{Int64,bicluster} )
#     global k_clust, iter, startTime
#     time_elapsed = (time() - startTime)/60
#     ## To get memory used -- DONE: add this to the stats_df; TODO: get this during MEME running (that's when max RAM is used) 
#     ## ps -U dreiss --no-headers -o rss -o comm | grep -E 'julia|meme|mast' | awk '{print $1}' | ( tr '\n' + ; echo 0 ) | bc
#     if OS_NAME != :Darwin
#         tmp = 0
#         try ## this command doesnt work on osiris although it works on fossil... weird
#             tmp = split( readall(`ps -U dreiss --no-headers -o rss -o comm` | `grep -E 'julia|meme|mast'` | `awk '{print $1}'`), '\n' )
#             tmp = sum([ parse_int(tmp[i]) for i=1:(length(tmp)-1) ])
#         catch
#             tmp = 0
#         end
#     else
#         tmp = split( readall(`ps -U dreiss -o rss -o comm` | `grep -E 'julia|meme|mast'` | `awk '{print $1}'`), '\n' )
#         tmp = sum([ parse_int(tmp[i]) for i=1:(length(tmp)-1) ])
#     end
#     (weight_r, weight_n, weight_m, weight_c, weight_v) = get_score_weights()
#     out_df = DataFrame( { "iter" => iter, "time" => time_elapsed, "mem_used" => tmp, 
#                          "r0" => weight_r, "n0" => weight_n, "m0" => weight_m, 
#                          "c0" => weight_c, "v0" => weight_v } )
#     tmp = float32([length(clusters[k].rows) for k=1:k_clust])
#     out_df["ROWS"] = nanmean(tmp)
#     println( "ROWS: ", out_df["ROWS"][1], " +/- ", nansd(tmp) )
#     tmp = float32([length(clusters[k].cols) for k=1:k_clust])
#     out_df["COLS"] = nanmean(tmp)
#     println( "COLS: ", out_df["COLS"][1], " +/- ", nansd(tmp) )
#     tmp = float32([clusters[k].resid for k=1:k_clust])
#     out_df["RESID"] = nanmean(tmp)
#     println( "MEAN RESID: ", out_df["RESID"][1], " +/- ", nansd(tmp) )
#     tmp = float32([clusters[k].dens_string for k=1:k_clust])
#     out_df["STRING_DENS"] = nanmean(tmp)
#     println( "MEAN STRING DENS: ", out_df["STRING_DENS"][1], " +/- ", nansd(tmp) )
#     tmp = float32([clusters[k].meanp_meme for k=1:k_clust])
#     out_df["MEME_PVAL"] = nanmean(tmp)
#     println( "MEAN MEME LOG10(P-VAL): ", out_df["MEME_PVAL"][1], " +/- ", nansd(tmp) )
#     rows = 0; for k in 1:k_clust rows = [rows, clusters[k].rows]; end
#     tmp = float32(collect(values(table(rows))))
#     out_df["CLUSTS_PER_ROW"] = nanmean(tmp)
#     println( "CLUSTS PER ROW: ", out_df["CLUSTS_PER_ROW"][1], " +/- ", nansd(tmp) )
#     cols = 0; for k in 1:k_clust cols = [cols, clusters[k].cols]; end
#     tmp = float32(collect(values(table(cols))))
#     out_df["CLUSTS_PER_COL"] = nanmean(tmp)
#     println( "CLUSTS PER COL: ", out_df["CLUSTS_PER_COL"][1], " +/- ", nansd(tmp) )
#     out_df
# end

# function clusters_to_dataFrame( clusters::Dict{Int64,bicluster} )
#     out = Array(DataFrame,length(clusters))
#     for k in 1:length(clusters)
#         b = clusters[k]
#         out_r = DataFrame( { "k" => k,
#                             "rows" => [join(rownames(ratios)[b.rows],',')],
#                             "resid" => b.resid,
#                             "dens_string" => b.dens_string,
#                             "meanp_meme" => b.meanp_meme,
#                             "cols" => [join(colnames(ratios)[b.cols],',')],
#                             "meme_out" => [join(b.meme_out,"<<<<>>>>")] 
#                             } )
#         out[k] = out_r
#     end
#     rbind(out)
# end
