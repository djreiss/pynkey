import sys
import copy

import numpy as np
from numpy import nan as NA
from numpy import random as rand
import pandas as pd

from colorama import Fore, Back, Style ## see https://pypi.python.org/pypi/colorama

from utils import slice_sampler
import scores
import funcs
import params

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
    meme_out = None ## string vector - meme output
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
        self.meme_out = ['']
        self.mast_out = pd.DataFrame()
        self.changed = np.ones(2, np.bool)

    def __repr__(self):
        return Fore.GREEN + 'Bicluster: ' + Fore.RESET + '%d' % self.k + '\n' + \
            '   resid: %f' % self.resid + '\n' + \
            '   meme-p: %f' % self.meanp_meme + '\n' + \
            '   string: %f' % self.dens_string + '\n' + \
            '   rows: %s' % str(self.rows) + '\n' + \
            '   cols: %s' % str(self.cols)
    
    def copy( self, deep=True ):
        out = copy.deepcopy(self) if deep else copy.copy(self)
        return out

    def volume( self ):
        return float(len(self.rows) * len(self.cols))

    def compute_residue( self, ratios ):
        rats = ratios.ix[ self.rows, self.cols ]
        ##if np.ndim( rats ) < 2 or np.size( rats, 0 ) <= 1 or np.size( rats, 1 ) <= 1 or \
        ##        np.mean( rats.isnull().values ) > 0.95:
        ##    warnings.warn( "COULD NOT COMPUTE RESIDUE FOR BICLUSTER " + self.k )
        ##    return 1.0
        return funcs.matrix_residue( rats )

    ## Note: how to profile this in ipython:
    ## %load_ext line_profiler
    ## %lprun -f funcs.matrix_residue Bicluster.bicluster.compute_residue_deltas(clust,ratios,all_genes)
    ##    or
    ## %lprun -f Bicluster.bicluster.compute_residue_deltas Bicluster.bicluster.compute_residue_deltas(clust,ratios,all_genes)
    def compute_residue_deltas( self, ratios, all_genes, actually_cols=False ):
        if not actually_cols:
            is_in = np.in1d( all_genes, self.rows )
            rows = self.rows
            rats = ratios.ix[ :, self.cols ]
            resid = funcs.matrix_residue( rats.ix[ rows, : ] )
        else:  ## all_genes is actually all_cols==ratios.columns.values
            is_in = np.in1d( all_genes, self.cols )
            rows = self.cols
            rats = ratios.ix[ self.rows, : ]
            resid = funcs.matrix_residue( rats.ix[ :, rows ] )
        all_resids = np.zeros( len( all_genes ), float )
        for i in xrange(len(all_genes)):
            r = all_genes[i]
            rows2 = rows[ rows != r ] if is_in[i] else np.append( rows, r )
            if not actually_cols:
                all_resids[i] = funcs.matrix_residue( rats.ix[ rows2, : ] )
            else:
                all_resids[i] = funcs.matrix_residue( rats.ix[ :, rows2 ] )
        return all_resids - resid

    def compute_var( self, ratios, var_add=0.1 ):
        rats = ratios.ix[ self.rows, self.cols ] ##.copy() ## get the bicluster's submatrix of the data
        return funcs.matrix_var( rats, var_add )

    def compute_var_deltas( self, ratios, all_genes, var_add=0.1, actually_cols=False ):
        if not actually_cols:
            is_in = np.in1d( all_genes, self.rows )
            rows = self.rows
            rats = ratios.ix[ :, self.cols ]
            resid = funcs.matrix_var( rats.ix[ rows, : ], var_add )
        else:  ## all_genes is actually all_cols==ratios.columns.values
            is_in = np.in1d( all_genes, self.cols )
            rows = self.cols
            rats = ratios.ix[ self.rows, : ]
            resid = funcs.matrix_var( rats.ix[ :, rows ], var_add )
        all_resids = np.zeros( len( all_genes ), float )
        for i in xrange(len(all_genes)):
            r = all_genes[i]
            rows2 = rows[ rows != r ] if is_in[i] else np.append( rows, r )
            if not actually_cols:
                all_resids[i] = funcs.matrix_var( rats.ix[ rows2, : ], var_add )
            else:
                all_resids[i] = funcs.matrix_var( rats.ix[ :, rows2 ], var_add )
        return all_resids - resid

    def compute_network_density( self, network ):
        return funcs.subnetwork_density( self.rows, network )

    def compute_network_density_deltas( self, network, all_genes ):
        is_in = np.in1d( all_genes, self.rows )
        rows = self.rows
        net1 = network.ix[ self.rows ]
        net1.set_index( ['protein2'], inplace=True )
        dens = funcs.subnetwork_density( rows, net1, already_subnetted=True )
        all_dens = np.zeros( len( all_genes ), float )
        for i in xrange(len(all_genes)):
            r = all_genes[i]
            rows2 = np.append( rows, r ) if is_in[i] else rows[ rows != r ]
            all_dens[i] = funcs.subnetwork_density( rows2, net1, already_subnetted=True )
        return all_dens - dens

    def compute_meme_pval( self ):
        if np.size(self.mast_out,0) <= 0:
            return NA
        df = self.mast_out.ix[ self.rows ] ## np.in1d( self.mast_out.Gene, self.rows ) ] ## make sure index of mast_out is Gene !!!
        mn = np.nanmean( log10( df['P-value'] ) )
        return mn

    ## Up-weight moves into clusters with <= 15 rows; down-weight moves out of clusters with >15 rows
    ## DONE: a little less heuristic?
    ## DONE? less "sharply-peaked" at 15? Perhaps make it "level out" between say 8 and 23
    ## DONE: plot this and see -- I think it DECREASES for really big volumes, want to fix this
    def get_volume_row_scores( self, all_genes ):
        thresh = params.avg_genes_per_cluster
        is_in = np.in1d( all_genes, self.rows )
        lr = len(self.rows)
        score_vr = np.array( [ (+1.0 if i else -1.0) * ( thresh - lr ) for i in is_in ] )
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
        ##score_vc = fill(NA, length(b.scores_c)) ## Do we really care how many cols there are?
        return score_vc

    ## Up-weight moves OUT if counts_g is HIGH, and moves IN if counts_g is LOW
    ## counts_g comes from counts_g = bicluster.get_all_cluster_row_counts( clusters, all_genes )
    def get_row_count_scores( self, all_genes, counts_g ):
        is_in = np.in1d( all_genes, self.rows )
        score_g = bicluster.get_cluster_row_count_scores( counts_g )
        score_g = np.array( [ (+score_g[all_genes[i]] if is_in[i] else -score_g[all_genes[i]]) for i in xrange(len(is_in)) ] )
        return score_g

    ## counts_g comes from counts_g = bicluster.get_all_cluster_row_counts( clusters, all_genes )
    ## TBD: check "changed" (rows=0, cols=1) and only recompute if True; then set "changed" to False. 
    def fill_all_scores(self, all_genes, ratios, string_net, counts_g, all_cols):
        if self.changed[0]:
            self.resid = self.compute_residue( ratios )
            self.dens_string = self.compute_network_density( string_net )
            self.scores_r = self.compute_residue_deltas( ratios, all_genes )
            self.scores_n = self.compute_network_density_deltas( string_net, all_genes )
            self.changed[0] = False
        if self.changed[1]:
            self.scores_c = self.compute_residue_deltas( ratios, all_cols, actually_cols=True )
            self.changed[1] = False
        return self

    ## counts_g comes from counts_g = bicluster.get_all_cluster_row_counts( clusters, all_genes )
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
        if length(self.rows) < min_rows or length(self.rows) > max_rows:
            nr = length(self.rows)
            warnings.warn( "RESEEDING BICLUSTER %d (%d)" % (self.k, nr) )

        ## DONE: add rows that are preferentially in few (or no) other clusters -- 
        ## sample from get_cluster_row_counts(clusters)
        ##clust.rows = unique([clust.rows, [Distributions.sample([1:nrow(ratios.x)]) for i=1:5]]) ## add 5 random rows
        counts_g = bicluster.get_all_cluster_row_counts( clusters, all_genes )
        g = np.array( counts_g.keys() )
        counts_g = np.array( counts_g.values() )
        if len(self.rows) < min_rows:
            counts_g = np.max(counts_g) + 0.01 - counts_g
            counts_g = counts_g / np.max(counts_g)
            self.rows = g[ unique( slice_sampler( counts_g, 5 ) ) ]
            ## add x/3 random cols
            self.cols = np.sort(np.unique(np.concatenate( (self.cols, \
                                                           ratios.columns.values[ rand.choice(ratios.shape[1], 
                                                           ratios.shape[1]/3-len(self.cols), replace=False) ] )), 0))
        elif len(self.rows) > max_rows:
            counts_g = ( counts_g + 0.01 ) / ( np.max(counts_g) + 0.01 )
            tmp_rows = g[np.in1d(g,self.rows)][ unique( slice_sampler( counts_g[np.in1d(g,self.rows)], 
                                                                        len(self.rows)-max_rows+1 ) ) ]
            self.rows = self.rows[ np.logical_not( in1d(self.rows, tmp_rows) ) ]
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

    ## THis is a static function so outside the definition of the bicluster
    @staticmethod
    def fill_all_cluster_scores(clusters, all_genes, ratios, string_net, all_conds):
        counts_g = bicluster.get_all_cluster_row_counts( clusters, all_genes )
        for k in clusters:
            print k
            clusters[k].fill_all_scores(all_genes, ratios, string_net, counts_g, ratios.columns.values)

    ## counts_g comes from counts_g = bicluster.get_all_cluster_row_counts( clusters, all_genes )
    @staticmethod
    def get_cluster_row_count_scores( counts_g ):
        thresh = params.avg_clusters_per_gene ##1.3 ## 2.0 ## 3.0 ## lower is better; coerce removing if gene is in more than 2 clusters
        return dict( { (i,j-thresh) for i,j in counts_g.items() } )

    ## TBD: write a generic update scores function that tests "update_func" (e.g. compute_resid) if a
    ##    each gene is added/removed from the cluster
    @staticmethod
    def get_update_scores( cluster, update_func, *args ):
        print args
        return update_func( cluster, *args )

##################################
## Note this works from ipython shell, but NOT as called in floc.get_floc_scores_all()
    ##import multiprocessing as mp
    import pathos.multiprocessing as mp

    @staticmethod
    def fill_all_cluster_scores_par(clusters, all_genes, ratios, string_net, all_conds, counts_g, threads=4):
        pool = mp.Pool(processes=threads)              # start 4 worker processes
        ## Need to send, e.g. a tuple (1, counts_g) if fill_all_scores_par() took multiple args
        clusters = pool.map(bicluster.fill_all_scores_par, clusters.keys() )
        pool.terminate()
        return clusters

    @staticmethod
    def fill_all_scores_par(k):
        print k
        clust = clusters[k]
        clust.fill_all_scores(all_genes, ratios, string_net, counts_g, ratios.columns.values)
        return clust
