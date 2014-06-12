import warnings
import datetime

import numpy as np
import pandas as pd

print 'importing funcs'

def copy_clusters( clusters, deep=True ):
    new_clusters = { k:clust.copy(deep=deep) for (k,clust) in clusters.items() }
    return new_clusters

def re_seed_all_clusters_if_necessary( clusters, ratios, all_genes, min_rows=3, max_rows=80 ):
    for k in clusters.keys():
        clusters[k] = clusters[k].re_seed_if_necessary( clusters, ratios, all_genes, min_rows, max_rows )
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

def print_cluster_stats( clusters, ratios, iter, startTime ):
    time_elapsed = datetime.datetime.now() - startTime ## seconds

    weight_r, weight_n, weight_m, weight_c, weight_v, weight_g = scores.get_score_weights( iter, ratios )

    out_df = DataFrame( { 'iter': iter, 'time': time_elapsed, ##'mem_used': tmp, 
                          'r0': weight_r, 'n0': weight_n, 'm0': weight_m, 
                          'c0': weight_c, 'v0': weight_v } )
    tmp = np.array([len(clusters[k].rows) for k in xrange(len(clusters))])
    out_df['ROWS'] = np.nanmean(tmp)
    print 'ROWS:', out_df.ROWS, ' +/- ', np.nanstd(tmp)
    tmp = np.array([len(clusters[k].cols) for k in xrange(len(clusters))])
    out_df['COLS'] = np.nanmean(tmp)
    print 'COLS:', out_df.COLS, ' +/- ', np.nanstd(tmp)
    tmp = np.array([clusters[k].resid for k in xrange(len(clusters))])
    out_df['RESID'] = np.nanmean(tmp)
    print 'MEAN RESID:', out_df.RESID, ' +/- ', np.nanstd(tmp)
    tmp = np.array([clusters[k].dens_string for k in xrange(len(clusters))])
    out_df['STRING_DENS'] = np.nanmean(tmp)
    print 'MEAN STRING DENS:', out_df.STRING_DENS, ' +/- ', nansd(tmp)
    tmp = np.array([clusters[k].meanp_meme for k in xrange(len(clusters))])
    out_df['MEME_PVAL'] = np.nanmean(tmp)
    print 'MEAN MEME LOG10(P-VAL):', out_df.MEME_PVAL, ' +/- ', nansd(tmp)
#     rows = 0; for k in 1:k_clust rows = [rows, clusters[k].rows]; end
#     tmp = np.array(collect(values(table(rows))))
#     out_df['CLUSTS_PER_ROW'] = nanmean(tmp)
#     print 'CLUSTS PER ROW:', out_df['CLUSTS_PER_ROW'][1], ' +/- ', nansd(tmp) )
#     cols = 0; for k in 1:k_clust cols = [cols, clusters[k].cols]; end
#     tmp = np.array(collect(values(table(cols))))
#     out_df['CLUSTS_PER_COL'] = nanmean(tmp)
#     print 'CLUSTS PER COL:', out_df['CLUSTS_PER_COL'][1], ' +/- ', nansd(tmp) )

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
