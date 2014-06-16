import warnings
import datetime

import numpy as np
import pandas as pd

import scores

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
    dens = float(np.sum(net2.weight)) / (float(len(rows))**2) ## Already symmetrized, need to decrease count by 1/2
    return np.log10( dens+1e-9 )

from collections import Counter
import resource

def print_cluster_stats( clusters, ratios, iter, startTime ):
    time_elapsed = datetime.datetime.now() - startTime ## seconds

    weight_r, weight_n, weight_m, weight_c, weight_v, weight_g = scores.get_score_weights( iter, ratios )

    ## To get memory used -- TODO: get this during MEME running (that's when max RAM is used) 
#     tmp = split( readall(`ps -U dreiss -o rss -o comm` | `grep -E 'julia|meme|mast'` | `awk '{print $1}'`), '\n' )
#     tmp = sum([ parse_int(tmp[i]) for i=1:(length(tmp)-1) ])
    ## see: https://stackoverflow.com/questions/938733/total-memory-used-by-python-process
    tmp = ( resource.getrusage(resource.RUSAGE_SELF).ru_maxrss + \
        resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss ) / 1024 ## in MB

    out_df = pd.DataFrame( { 'iter': [iter], 'time': [time_elapsed], 'mem_used': tmp, 
                          'r0': [weight_r], 'n0': [weight_n], 'm0': [weight_m], 
                          'c0': [weight_c], 'v0': [weight_v] } )
    tmp = np.array([len(clusters[k].rows) for k in xrange(len(clusters))])
    out_df['ROWS'] = np.nanmean(tmp)
    print 'ROWS:', out_df.ROWS[0], ' +/- ', np.nanstd(tmp)
    tmp = np.array([len(clusters[k].cols) for k in xrange(len(clusters))])
    out_df['COLS'] = np.nanmean(tmp)
    print 'COLS:', out_df.COLS[0], ' +/- ', np.nanstd(tmp)
    tmp = np.array([clusters[k].resid for k in xrange(len(clusters))])
    out_df['RESID'] = np.nanmean(tmp)
    print 'MEAN RESID:', out_df.RESID[0], ' +/- ', np.nanstd(tmp)
    tmp = np.array([clusters[k].dens_string for k in xrange(len(clusters))])
    out_df['STRING_DENS'] = np.nanmean(tmp)
    print 'MEAN STRING DENS:', out_df.STRING_DENS[0], ' +/- ', np.nanstd(tmp)
    tmp = np.array([clusters[k].meanp_meme for k in xrange(len(clusters))])
    tmp[ np.isinf(tmp) ] = np.nan
    out_df['MEME_PVAL'] = np.nanmean(tmp)
    print 'MEAN MEME LOG10(P-VAL):', out_df.MEME_PVAL[0], ' +/- ', np.nanstd(tmp)
    rows = clusters[0].rows
    for k in range(1,len(clusters)):
        rows = np.append(rows, clusters[k].rows)
    c = Counter(rows)
    tmp = np.array(c.values())
    out_df['CLUSTS_PER_ROW'] = np.nanmean(tmp)
    print 'CLUSTS PER ROW:', out_df.CLUSTS_PER_ROW[0], ' +/- ', np.nanstd(tmp)
    cols = clusters[0].cols
    for k in range(1,len(clusters)):
        cols = np.append(cols, clusters[k].cols)
    c = Counter(cols)
    tmp = np.array(c.values())
    out_df['CLUSTS_PER_COL'] = np.nanmean(tmp)
    print 'CLUSTS PER COL:', out_df.CLUSTS_PER_COL[0], ' +/- ', np.nanstd(tmp)
    return out_df

def clusters_to_dataFrame( clusters ):
    out = {}
    for k in xrange(1,len(clusters)):
        b = clusters[k]
        out_r = pd.DataFrame( { 'k': [k],
                                'rows': ','.join(b.rows),
                                'resid': [b.resid],
                                'dens_string': [b.dens_string],
                                'meanp_meme': [b.meanp_meme],
                                'cols': ','.join(b.cols),
                                'meme_out': '<<<<>>>>'.join('\n'.split(b.meme_out)) } )
        out[k] = out_r
    out = pd.concat( out.values() )
    return out

import re

## Find the flagellar cluster, whew!!!
def clusters_w_func( func, clusters, anno, n_best=1 ):
    reg = re.compile(func)
    inds = np.nonzero( np.array( [ len( reg.findall( str(i) ) ) for i in anno.desc.values ] ) )
    genes = anno.index.values[ inds ]
    nhits = np.array( [ np.sum( np.in1d( genes, b.rows ) ) for b in clusters.values() ] )
    ords = np.argsort( nhits )
    kInds = ords[ len(ords)-n_best: ]
    
    for kInd in kInds:
        genes = clusters[kInd].rows
        ##print genes ## print the genes
        genes = genes[ np.in1d(genes, anno.index.values) ]
        print kInd, nhits[kInd], '\n', anno.ix[ np.in1d(anno.index.values,genes), ['desc'] ]
    return kInds

