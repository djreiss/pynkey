import copy
import warnings

import numpy as np
from numpy import nan as NA
from numpy import random as rand
import pandas as pd

from colorama import Fore, Back, Style ## see https://pypi.python.org/pypi/colorama

import scores
import funcs
import params
import utils as ut ##slice_sampler, do_something_par
import sequence as seq
import meme

print 'importing Bicluster'

class bicluster:
    k = None    ## cluster index
    rows = None ## vector of gene names
    cols = None ## vector of cond names
    var = None  ## float - variance??
    resid = None ## float - residual
    dens_string = None ## float - string network density
    meanp_meme = None ## float - mean meme p-value
    scores_r = None ## vector of row scores
    scores_c = None ## vector of col scores
    scores_n = None ## vector of network scores
    scores_m = None ## vector of motif scores
    meme_out = None ## big string (containing \n's) - meme output
    mast_out = None ## DataFrame - parsed meme output
    changed = None

    def __init__( self, k, rows, cols=None, ratios=None ):
        self.k = k
        self.rows = np.unique(rows)
        if cols is not None:
            self.cols = np.unique(cols)
        elif ratios is not None:
            self.cols = np.sort( ratios.columns.values[ rand.choice(ratios.shape[1], ratios.shape[1]/2, replace=False) ] )

        max_float = float('inf') ##sys.float_info.max
        self.var = max_float
        self.resid = max_float
        self.dens_string = max_float
        self.meanp_meme = max_float
        self.scores_r = np.array([])
        self.scores_c = np.array([])
        self.scores_n = np.array([])
        self.scores_m = np.array([])
        self.meme_out = ''
        self.mast_out = pd.DataFrame()
        self.changed = np.ones(2, np.bool)

    def __repr__(self):
        return Fore.GREEN + 'Bicluster: ' + Fore.RESET + '%d' % self.k + '\n' + \
            '   resid: %f' % self.resid + '\n' + \
            '   meme-p: %f' % self.meanp_meme + '\n' + \
            '   string: %f' % self.dens_string + '\n' + \
            '   rows: %s' % str(self.rows) + '\n' + \
            '   cols: %s' % str(self.cols) + '\n' + \
            '   changed: %s' % str(self.changed)
    
    def copy( self, deep=True ):
        out = copy.deepcopy(self) if deep else copy.copy(self)
        return out

    def volume( self ):
        return float(len(self.rows) * len(self.cols))

    def compute_residue( self, ratios ):
        if len(self.rows) <= 0:
            return NA
        rats = ratios.ix[ self.rows, self.cols ]
        ##if np.ndim( rats ) < 2 or np.size( rats, 0 ) <= 1 or np.size( rats, 1 ) <= 1 or \
        ##        np.mean( rats.isnull().values ) > 0.95:
        ##    warnings.warn( "COULD NOT COMPUTE RESIDUE FOR BICLUSTER " + self.k )
        ##    return 1.0
        return funcs.matrix_residue( rats.values )

    ## Note: how to profile this in ipython:
    ## %load_ext line_profiler
    ## %lprun -f funcs.matrix_residue Bicluster.bicluster.compute_residue_deltas(clust,ratios,all_genes)
    ##    or
    ## %lprun -f Bicluster.bicluster.compute_residue_deltas Bicluster.bicluster.compute_residue_deltas(clust,ratios,all_genes)
    def compute_residue_deltas( self, ratios, all_genes, actually_cols=False ):
        if not actually_cols:
            if len(self.rows) <= 0:
                return np.repeat( -0.1, len(all_genes) )
            in_rats = ratios
            rows = self.rows
            cols = self.cols
            is_inRats = np.in1d( all_genes, ratios.index.values )
        else:  ## all_genes needs to be all_cols==ratios.columns.values
            if len(self.cols) <= 0:
                return np.repeat( -0.1, len(all_genes) )
            in_rats = ratios.T
            rows = self.cols
            cols = self.rows
            is_inRats = np.in1d( all_genes, ratios.columns.values )
        is_in = np.in1d( all_genes, rows )
        rats = in_rats.ix[ :, cols ]
        resid = funcs.matrix_residue( rats.ix[ rows ].values )
        all_resids = np.zeros( len( all_genes ), float )
        for i in xrange(len(all_genes)):
            if not is_inRats[i]:
                all_resids[i] = resid
                continue
            r = all_genes[i]
            rows2 = rows[ rows != r ] if is_in[i] else np.append( rows, r )
            all_resids[i] = funcs.matrix_residue( rats.ix[ rows2 ].values )
        return all_resids - resid

    ## 143 ms vs. 355 ms for old way
    ## Faster way seems to be to take bicluster submatrix, then append single row for 
    ##   each test (rather than old way which is to subset entire matrix each iter)
    def compute_residue_deltas_FASTER( self, ratios, all_genes, actually_cols=False ):
        ## if actually_cols is False, then all_genes should be ratios.colums.values
        if not actually_cols:
            if len(self.rows) <= 0:
                return np.repeat( -0.1, len(all_genes) )
            in_rats = ratios
            rows = self.rows
            cols = self.cols
            is_inRats = np.in1d( all_genes, ratios.index.values )
        else:  ## all_genes needs to be all_cols==ratios.columns.values
            if len(self.cols) <= 0:
                return np.repeat( -0.1, len(all_genes) )
            in_rats = ratios.T
            rows = self.cols
            cols = self.rows
            is_inRats = np.in1d( all_genes, ratios.columns.values )

        is_in = np.in1d( all_genes, rows )
        rats = in_rats.ix[ :, cols ]
        rrats = rats.ix[ rows ]
        resid = funcs.matrix_residue( rrats.values )
        all_resids = np.zeros( len( all_genes ), float )
        for i in xrange(len(all_genes)):
            if not is_inRats[i]:
                all_resids[i] = resid
                continue
            r = all_genes[i]
            if is_in[i]:
                rows2 = rows[ rows != r ]
                tmp_rats = rrats.ix[ rows2 ].values
            else:
                tmp_rats = np.append( rrats.values, [rats.ix[ r ].values], axis=0 )
            all_resids[i] = funcs.matrix_residue( tmp_rats )
        return all_resids - resid

    def compute_var( self, ratios, var_add=0.1 ):
        rats = ratios.ix[ self.rows, self.cols ] ##.copy() ## get the bicluster's submatrix of the data
        return funcs.matrix_var( rats, var_add )

    ## TBD: re-design this function as is done by compute_residue_deltas
    def compute_var_deltas( self, ratios, all_genes, var_add=0.1, actually_cols=False ):
        if not actually_cols:
            is_in = np.in1d( all_genes, self.rows )
            rows = self.rows
            rats = ratios.ix[ :, self.cols ]
            resid = funcs.matrix_var( rats.ix[ self.rows ], var_add )
        else:  ## all_genes is actually all_cols==ratios.columns.values
            is_in = np.in1d( all_genes, self.cols )
            rows = self.cols
            rats = ratios.ix[ self.rows ]
            resid = funcs.matrix_var( rats.ix[ :, self.rows ], var_add )
        all_resids = np.zeros( len( all_genes ), float )
        for i in xrange(len(all_genes)):
            r = all_genes[i]
            rows2 = rows[ rows != r ] if is_in[i] else np.append( rows, r )
            if not actually_cols:
                all_resids[i] = funcs.matrix_var( rats.ix[ rows2 ], var_add )
            else:
                all_resids[i] = funcs.matrix_var( rats.ix[ :, rows2 ], var_add )
        return all_resids - resid

    def compute_network_density( self, network ):
        if len(self.rows) <= 0:
            return NA
        return funcs.subnetwork_density( self.rows, network )

    def compute_network_density_deltas( self, network, all_genes ):
        if len(self.rows) <= 0:
            return np.repeat( -0.1, len(all_genes) )
        is_in = np.in1d( all_genes, self.rows )
        rows = self.rows
    #     dens = funcs.subnetwork_density( rows, network )
        net1 = network[ network[[0,1]].isin(rows).any(1) ] ## subnetwork where protein1 OR protein2 are in rows
        if net1.shape[0] <= 0:
            return np.repeat( -0.1, len(all_genes) )
        dens = funcs.subnetwork_density( rows, net1 ) ## pass subnetwork to this func for speed
        all_dens = np.zeros( len( all_genes ), float )
        for i in xrange(len(all_genes)):
            r = all_genes[i]
            rows2 = rows[ rows != r ] if is_in[i] else np.append( rows, r )
            all_dens[i] = funcs.subnetwork_density( rows2, net1 )
        return all_dens - dens

    def get_sequences( self, anno, genome_seqs, op_table, distance, do_filter=True):
        seqs = seq.get_sequences(self.rows, anno, genome_seqs, op_table, True, distance, False) 
        if do_filter:
            seqs = seq.filter_sequences(seqs, distance)
        return seqs

    ## put it here instead of in meme.py, b/c its better here
    def re_meme( self, distance, allSeqs_fname, anno, genome_seqs, op_table, motif_width_range, n_motifs=2 ): 
        seqs = self.get_sequences( anno, genome_seqs, op_table, distance, do_filter=True )
        k, meme_out, mast_out = meme.re_meme_bicluster( self.k, seqs, n_motifs, allSeqs_fname, 
                                                        motif_width_range, False )
        if k != self.k:
            warnings.warn( 'Uh oh! %d' %k )
        else:
            self.meme_out = meme_out
            if mast_out.shape[0] > 0:
                mast_out = mast_out.set_index( 'Gene' )
            self.mast_out = mast_out
            self.changed[0] = True
        return self

    def compute_meme_pval( self ):
        if np.size(self.mast_out,0) <= 0 or len(self.rows) <= 0:
            return NA
        df = self.mast_out.ix[ self.rows ] ## make sure index of mast_out is Gene !!! Done - in bicluster.re_meme()
        mn = np.nanmean( np.log10( df['P-value'] ) )
        return mn

    def compute_meme_pval_deltas( self, all_genes ):
        if len(self.rows) <= 0:
            return np.repeat( -0.1, len(all_genes) )
        is_in = np.in1d( all_genes, self.rows )
        pval = self.compute_meme_pval()
        all_pvals = np.zeros( len( all_genes ), float )
        for i in xrange(len(all_genes)):
            r = all_genes[i]
            rows2 = self.rows[ self.rows != r ] if is_in[i] else np.append( self.rows, r ) 
            df = self.mast_out.ix[ rows2 ] ## make sure index of mast_out is Gene !!! Done - in bicluster.re_meme()
            all_pvals[i] = np.nanmean( np.log10( df['P-value'] ) )
        return all_pvals - pval

    ## Up-weight moves into clusters with <= 15 rows; down-weight moves out of clusters with >15 rows
    ## DONE: a little less heuristic?
    ## DONE? less "sharply-peaked" at 15? Perhaps make it "level out" between say 8 and 23
    ## DONE: plot this and see -- I think it DECREASES for really big volumes, want to fix this
    def get_volume_row_scores( self, all_genes ):
        thresh = params.avg_genes_per_cluster
        is_in = np.in1d( all_genes, self.rows )
        lr = len(self.rows)
        score_vr = np.array( [ (+1.0 if i else -1.0) * ( thresh - lr ) for i in is_in ] )
        score_vr = score_vr**3.0
        if lr >= thresh - 7 and lr <= thresh + 7:
            score_vr /= 5.0
        elif lr >= thresh - 12 and lr <= thresh + 12:
            score_vr /= 2.0
        return score_vr

    ## Just up-weight moves that add columns (to prevent shrinkage)
    ## NOTE that as is, the weighting decreases for really high #cols... is this what we want?
    def get_volume_col_scores( self, all_cols ):
        is_in = np.in1d( all_cols, self.cols )
        lc = len(self.cols)
        score_vc = np.array( [ +1.0/lc if i else -1.0/lc for i in is_in ] )
        ##score_vc = fill(NA, len(b.scores_c)) ## Do we really care how many cols there are?
        return score_vc

    ## Up-weight moves OUT if counts_g is HIGH, and moves IN if counts_g is LOW
    ## counts_g comes from counts_g = bicluster.get_all_cluster_row_counts( clusters, all_genes )
    def get_row_count_scores( self, all_genes, counts_g ):
        is_in = np.in1d( all_genes, self.rows )
        score_g = bicluster.get_cluster_row_count_scores( counts_g )
        score_g2 = np.array( [ (-score_g[all_genes[i]] if is_in[i] else +score_g[all_genes[i]]) \
                                  for i in xrange(len(is_in)) ] )
        return score_g2

    ## counts_g comes from counts_g = bicluster.get_all_cluster_row_counts( clusters, all_genes )
    ## TBD: check "changed" (rows=0, cols=1) and only recompute if True; then set "changed" to False. 
    ## Note residual stuff needs to be computed if *either* rows *or* columns have changed.
    def fill_all_scores(self, iter, all_genes, ratios, string_net, counts_g, all_cols, force=False):
        weight_r, weight_n, weight_m, weight_c, weight_v, weight_g = scores.get_score_weights( iter, ratios )
        if not force and np.all(np.invert(self.changed)):
            return self
        print self.k, self.changed
        if weight_r > 0:
            self.resid = self.compute_residue( ratios )  ## do this if *either* rows or cols changed
            self.scores_r = self.compute_residue_deltas( ratios, all_genes )
            self.scores_c = self.compute_residue_deltas( ratios, all_cols, actually_cols=True )
        if self.changed[0] or force:
            if weight_m > 0 and self.meme_out != '':
                self.meme_pval = self.compute_meme_pval()
                self.scores_m = self.compute_meme_pval_deltas( all_genes )
            if abs(weight_n) > 0:
                self.dens_string = self.compute_network_density( string_net )
                self.scores_n = self.compute_network_density_deltas( string_net, all_genes )
        self.changed[0] = self.changed[1] = False
        return self

    ## counts_g comes from counts_g = bicluster.get_all_cluster_row_counts( clusters, all_genes )
    ## Get DataFrame formatted with all scores, for floc.update()
    def get_floc_scoresDF_rows(self, iter, all_genes, ratios, counts_g):
        is_in = np.in1d( all_genes, self.rows )
        ## only use them if their weights are > 0 and they have been filled via fill_all_scores()
        weight_r, weight_n, weight_m, weight_c, weight_v, weight_g = scores.get_score_weights( iter, ratios )
        NAs =  np.repeat(NA, len(self.scores_r)) if (weight_r <= 0 or weight_n <= 0 or weight_m <= 0) else NA
        score_r = self.scores_r if weight_r > 0 and len(self.scores_r)>0 else NAs
        score_n = self.scores_n if abs(weight_n) > 0 and len(self.scores_n)>0 else NAs
        score_m = self.scores_m if weight_m > 0 and len(self.scores_m)>0 else NAs
        score_vr = self.get_volume_row_scores( all_genes )
        score_g = self.get_row_count_scores( all_genes, counts_g )
        out = pd.DataFrame( { 'row_col':all_genes,
                              'is_in':is_in,
                              'is_row_col':np.repeat('r', len(self.scores_r)), ## CANT: move this outside the loop
                              'k':np.repeat(self.k, len(self.scores_r)),
                              'score_r':score_r,
                              'score_n':score_n,
                              'score_m':score_m,
                              'score_v':score_vr,
                              'score_g':score_g } )
        return out

    def get_floc_scoresDF_cols(self, iter, ratios):
        all_cols = ratios.columns.values
        is_in = np.in1d( all_cols, self.cols )
        (weight_r, weight_n, weight_m, weight_c, weight_v, weight_g) = scores.get_score_weights( iter, ratios )
        NAs = np.repeat(NA, len(self.scores_c)) ## if weight_c <= 0 else NA
        score_c = self.scores_c if weight_c > 0 else NAs
        score_vc = self.get_volume_col_scores( all_cols )
        out = pd.DataFrame( { 'row_col':all_cols,
                              'is_in':is_in,
                              'is_row_col':np.repeat('c', len(self.scores_c)), ## CANT: move this outside the loop
                              'k':np.repeat(self.k, len(self.scores_c)),
                              'score_r':score_c,
                              'score_n':NAs, ## CANT: move this outside the loop
                              'score_m':NAs,  ## CANT: move this outside the loop
                              'score_v':score_vc,
                              'score_g':NAs } )
        return out

    def re_seed_if_necessary( self, clusters, ratios, all_genes, min_rows=3, max_rows=80 ):
        ## Need to fill in "empty" clusters, mostly because var is non-defined, Also because meme-ing only works 
        ## with >2 sequences. 
        ## Squash clusters that get way too big (DONE: just remove some to bring it down to max_rows)
        if len(self.rows) < min_rows or len(self.rows) > max_rows:
            nr = len(self.rows)
            warnings.warn( "RESEEDING BICLUSTER %d (%d)" % (self.k, nr) )

        ## DONE: add rows that are preferentially in few (or no) other clusters -- 
        ## sample from get_cluster_row_counts(clusters)
        ##clust.rows = np.unique([clust.rows, [Distributions.sample([1:nrow(ratios.x)]) for i=1:5]]) ## add 5 random rows
        counts_g = bicluster.get_all_cluster_row_counts( clusters, all_genes )
        g = np.array( counts_g.keys() )
        counts_g = np.array( counts_g.values() )
        if len(self.rows) < min_rows:
            counts_g = np.max(counts_g) + 0.01 - counts_g
            counts_g = counts_g / np.max(counts_g)
            self.rows = g[ np.unique( ut.slice_sampler( counts_g, 5 ) ) ]
            ## add x/3 random cols
            self.cols = np.sort(np.unique(np.concatenate( (self.cols, \
                                                           ratios.columns.values[ rand.choice(ratios.shape[1], \
                                                           ratios.shape[1]/3-len(self.cols), replace=False) ] )), 0))
            self.changed[0] = self.changed[1] = True
        elif len(self.rows) > max_rows:
            counts_g = ( counts_g + 0.01 ) / ( np.max(counts_g) + 0.01 )
            tmp_rows = g[np.in1d(g,self.rows)][ np.unique( ut.slice_sampler( counts_g[np.in1d(g,self.rows)], 
                                                                        len(self.rows)-max_rows+1 ) ) ]
            self.rows = self.rows[ np.logical_not( np.in1d(self.rows, tmp_rows) ) ]
            self.changed[ 0 ] = True
        return self

    @staticmethod
    ## number of clusters each gene is in - need to compute only once over all clusters
    ## all_genes provides reference of all possible gene names. May want to use ratios.index.values for this
    ## TBD NOW: instead of returning a dict, return a pandas Series. No, dict is more than 10x faster
    def get_all_cluster_row_counts( clusters, all_genes ):
        ##counts = np.array( [len(clust.rows) for clust in clusters.values()] )
        ##return np.bincount( counts )
        ##d = pd.Series(np.zeros(len(all_genes), int), all_genes) 
        d = dict( zip(list(all_genes), list(np.zeros(len(all_genes), int))) )
        for cc in clusters.values():
            for r in cc.rows:
                d[r] += 1
        return d

    ## counts_g comes from counts_g = bicluster.get_all_cluster_row_counts( clusters, all_genes )
    @staticmethod
    def get_cluster_row_count_scores( counts_g ):
        thresh = params.avg_clusters_per_gene ##1.3 ## 2.0 ## 3.0 ## lower is better; coerce removing if gene is in more than 2 clusters
        pow_to_use = 1.0 ## 3.0
        return dict( { (i,(j-thresh)**pow_to_use) for i,j in counts_g.items() } )

    ## TBD: write a generic update scores function that tests "update_func" (e.g. compute_resid) if a
    ##    each gene is added/removed from the cluster
    @staticmethod
    def get_update_scores( cluster, update_func, *args ):
        print args
        return update_func( cluster, *args )

    ## THis is a static function so outside the definition of the bicluster
    # @staticmethod
    # def fill_all_cluster_scores(clusters, iter, all_genes, ratios, string_net, all_conds):
    #     counts_g = bicluster.get_all_cluster_row_counts( clusters, all_genes )
    #     for k in clusters:
    #         print 'FILL: %s' % k
    #         clusters[k].fill_all_scores(iter, all_genes, ratios, string_net, counts_g, ratios.columns.values)
    #     return clusters

    ## THis is a static function so outside the definition of the bicluster
    # @staticmethod
    # def re_meme_all_clusters(clusters, distance_search, allSeqs_fname, anno, genome_seqs, op_table, 
    #                          motif_width_range, n_motifs=2 ): 
    #     for k in clusters:
    #         print 'MEME: %s' % k
    #         clusters[k].re_meme( distance_search, allSeqs_fname, anno, genome_seqs, op_table, 
    #                              motif_width_range, n_motifs )
    #     return clusters

##################################
## Note this works from ipython shell, but NOT as called in floc.get_floc_scores_all()
    # @staticmethod
    # def fill_all_cluster_scores_par( threads=4 ): ##clusters, all_genes, ratios, string_net, all_conds, counts_g, threads=4):
    #     global clusters, all_genes, ratios, string_net, counts_g
    #     pool = mp.Pool(processes=threads)              # start 4 worker processes
    #     ## Need to send, e.g. a tuple (1, counts_g) if fill_all_scores_par() took multiple args
    #     clusters = pool.map(bicluster.fill_all_scores_par, clusters.keys())
    #     pool.terminate()
    #     clusters = {clusters[i].k: clusters[i] for i in xrange(len(clusters))}  # convert back to map
    #     return clusters

    # @staticmethod
    # def fill_all_scores_par( k ):
    #     global clusters, all_genes, ratios, string_net, counts_g, ratios
    #     print k
    #     clust = clusters[k]
    #     clust.fill_all_scores(all_genes, ratios, string_net, counts_g, ratios.columns.values)
    #     return clust

import globals as glb

def fill_all_cluster_scores_par( clusters, threads=None, force=False ):
    ## Clusters per gene counts - precompute:
    glb.counts_g = bicluster.get_all_cluster_row_counts( clusters, glb.all_genes )
    if force:
        for c in clusters.values():
            c.changed[0] = c.changed[1] = True
    clusters = ut.do_something_par( clusters.values(), fill_all_scores_par, threads=threads )
    clusters = {clusters[i].k: clusters[i] for i in xrange(len(clusters))}  # convert back to map
    return clusters

## Keep everything global so it doesn't need to be sent to each child.
def fill_all_scores_par( clust ):
    weight_r, weight_n, weight_m, weight_c, weight_v, weight_g = scores.get_score_weights( glb.iter, glb.ratios )
    was_changed = clust.changed[0]
    if weight_m > 0 and was_changed:
        n_mots = meme.get_n_motifs( glb.iter, params.n_iters )
        clust.re_meme( params.distance_search, glb.allSeqs_fname, glb.anno, glb.genome_seqs, 
                       glb.op_table, params.motif_width_range, n_motifs=n_mots )
    clust.fill_all_scores(glb.iter, glb.all_genes, glb.ratios, glb.string_net, glb.counts_g, 
                          glb.ratios.columns.values)
    clust.changed[0] = clust.changed[1] = False
    return clust

## Now this only sends over k to the children; clusters and genome_seqs, etc. are all global in childrens' namespace
# def re_meme_all_clusters_par( clusters, threads=None ):
#     meme_outs = do_something_par( clusters.keys(), re_meme_par, threads=threads )
#     ##meme_outs = { i[0]: (i[1], i[2]) for i in meme_outs }
#     for i in meme_outs:
#         clusters[ i[0] ].meme_out = i[1]
#         clusters[ i[0] ].mast_out = i[2]
#     return clusters

# def re_meme_par( clustK ):
#     n_motifs = meme.get_n_motifs( glb.iter, params.n_iters )
#     clust = glb.clusters[ clustK ]
#     if not clust.changed[0]:
#         return (clustK, clust.meme_out, clust.mast_out 
#     seqs = clust.get_sequences( glb.anno, glb.genome_seqs, glb.op_table, params.distance_search, 
#                                 do_filter=True )
#     k, meme_out, mast_out = meme.re_meme_bicluster( clustK, seqs, n_motifs, glb.allSeqs_fname, 
#                                                     params.motif_width_range, verbose=False )
#     return (k, meme_out, mast_out)
