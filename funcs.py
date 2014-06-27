import numpy as np
import pandas as pd

import scores
import weaved

print 'importing funcs'

def copy_clusters( clusters, deep=True ):
    new_clusters = { k:clust.copy(deep=deep) for (k,clust) in clusters.items() }
    return new_clusters

def re_seed_all_clusters_if_necessary( clusters, ratios, all_genes, min_rows=3, max_rows=80 ):
    for k in clusters.keys():
        clusters[k] = clusters[k].re_seed_if_necessary( clusters, ratios, all_genes, min_rows, max_rows )
    return clusters

## rats needs to be a numpy array (not dataframe - so use df.values if necessary)
def matrix_residue( rats, weaveIt=True ):
    if np.ndim( rats ) < 2 or np.size( rats, 0 ) <= 1 or np.size( rats, 1 ) <= 1: ## or \
            ##np.mean( rats.isnull().values ) > 0.95:
        import logging
        logging.warning( "COULD NOT COMPUTE RESIDUE" )
        return 1.0

    ##rats = rats.values ## do it all in numpy - faster

    if not weaveIt:
        d_rows = np.nanmean(rats, 1) ##rats.mean(1)
        d_cols = np.nanmean(rats, 0) ##rats.mean(0)
        d_all = np.nanmean(d_rows) ##d_rows.mean() ## pandas default is to ignore nan's
        rats = rats + d_all - np.add.outer( d_rows, d_cols )
        ## another way of doing this, but about the same speed:
        ##rats = (rats+d_all-d_cols).T-d_rows
        average_r = np.nanmean( np.abs(rats) )

    average_r = weaved.fast_resid(rats)
    return average_r

def matrix_var( rats, var_add=0.1 ):
    rats = rats.values
    mn = np.nanmean(rats, 0) ## subtract bicluster mean profile
    return np.nanvar(rats-mn) / (np.nanvar(mn) + var_add)

# def subnetwork_density_OLD( rows, network ):
#     net1 = network.ix[rows]
#     net1.set_index( ['protein2'], inplace=True ) ##[ np.in1d( net1.protein2, rows ) ]
#     net2 = net1.ix[ rows ]
#     dens = float(np.sum(net2.weight)) / (float(len(rows))**2) ## Already symmetrized, need to decrease count by 1/2
#     return np.log10( dens+1e-9 )

## Use edge density score: sum of edge weights (=number of edges for all weights =1) / number of nodes^2
## This is the faster version that uses SubDataFrames
## TODO: use numexpr to speed up and avoid temporary array creation?
## This is slower but for accuracy when passed a subnetwork it can be faster
def subnetwork_density( rows, network ):
    net1 = network[ network[[0,1]].isin(rows).all(1) ]
    dens = float(np.sum(net1.weight)) / (float(len(rows))**2) ## Already symmetrized, need to decrease count by 1/2
    return np.log10( dens+1e-9 )

from collections import Counter
import resource

def print_cluster_stats( clusters, ratios, iter, startTime, n_tries=0, n_improvements=0, 
                         changed_rows=0, changed_cols=0 ):
    weight_r, weight_n, weight_m, weight_c, weight_v, weight_g = scores.get_score_weights( iter, ratios )
    print 'ITER:', iter
    print 'r0: %.3f; n0: %.3f; m0: %.3f; c0: %.3f; v0: %.3f, g0: %.3f' % ( weight_r, weight_n, weight_m, 
                                                                           weight_c, weight_v, weight_g )
    print 'N_MOVES:', n_tries
    print 'N_IMPROVEMENTS:', n_improvements

    import datetime
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
                          'c0': [weight_c], 'v0': [weight_v], 'g0': [weight_g] } )
    tmp = np.array([len(clusters[k].rows) for k in xrange(len(clusters))])
    out_df['ROWS'] = np.nanmean(tmp[tmp>1])
    print 'ROWS:', out_df.ROWS[0], ' +/- ', np.nanstd(tmp[tmp>1])
    tmp = np.array([len(clusters[k].cols) for k in xrange(len(clusters))])
    out_df['COLS'] = np.nanmean(tmp[tmp>1])
    print 'COLS:', out_df.COLS[0], ' +/- ', np.nanstd(tmp[tmp>1])
    tmp = np.array([clusters[k].resid for k in xrange(len(clusters))])
    out_df['RESID'] = np.nanmean(tmp[tmp>0.001])
    print 'MEAN RESID:', out_df.RESID[0], ' +/- ', np.nanstd(tmp[tmp>0.001])
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

    out_df['N_MOVES'] = n_tries
    out_df['N_IMPROVEMENTS'] = n_improvements
    out_df['N_CLUSTS_CHANGED_ROWS'] = changed_rows
    out_df['N_CLUSTS_CHANGED_COLS'] = changed_cols
    print 'N_CLUSTS_CHANGED (ROWS): ', out_df.N_CLUSTS_CHANGED_ROWS[0]
    print 'N_CLUSTS_CHANGED (COLS): ', out_df.N_CLUSTS_CHANGED_COLS[0]
    return out_df

## save out everything in params.* and globals.*
def checkpoint( fname, verbose=False ):
    import params,globals ## for loading/saving environment
    import cPickle as pickle
    import gzip
    d = dict() ## first build up a big dict with all items to be saved
    for k,v in params.__dict__.items(): ## save everything from params
        if k.startswith('__'): continue;
        if ( type(v) == type(pickle) ): continue; ## module type
        if verbose: 
            print 'Saving:', k, type(k)
        d[ 'params.%s'%k ] = v
    for k,v in globals.__dict__.items(): ## save everything from params
        if k.startswith('__'): continue;
        if ( type(v) == type(pickle) ): continue; ## module type
        if verbose:
            print 'Saving:', k, type(v)
        d[ 'globals.%s'%k ] = v
    f = gzip.open( fname, 'wb' )
    pickle.dump( d, f, pickle.HIGHEST_PROTOCOL ) ## dump the dict
    f.close()

def load_checkpoint( fname, verbose=False ): ## use 'exec' to load the values
    import params, globals ## for loading/saving environment ... need to allow for loading globals without running init
    ## trick from http://lucumr.pocoo.org/2011/2/1/exec-in-python/
    import cPickle as pickle
    import gzip
    f = gzip.open( fname, 'r' )
    dd = pickle.load( f ) ## read the dict
    f.close()
    for k,v in dd.items():
        if verbose:
            print 'Loading: %s' % k
        code = compile('%s = v' % k, '<string>', 'exec')
        exec code

def clusters_to_dataFrame( clusters ):
    out = {}
    for k in xrange(0,len(clusters)):
        b = clusters[k]
        out_r = pd.DataFrame( { 'k': [k+1],
                                'rows': ','.join(b.rows),
                                'resid': [b.resid],
                                'dens_string': [b.dens_string],
                                'meanp_meme': [b.meanp_meme],
                                'cols': ','.join(b.cols),
                                'meme_out': '' if b.meme_out == '' else '<<<<>>>>'.join(b.meme_out.split('\n')) } )
        out[k] = out_r
    out = pd.concat( out.values() )
    return out

import re

## Find the flagellar cluster, whew!!!
def clusters_w_func( func, clusters, anno, n_best=1 ):
    reg = re.compile(func, flags=re.IGNORECASE)
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
