import sys
import copy

from numpy import nan as NA
import numpy as np
from numpy import random as rand
import pandas as pd

from colorama import Fore, Back, Style ## see https://pypi.python.org/pypi/colorama

from slice_sampler import slice_sampler

import funcs

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
    changed = False 

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

    def residue( self, ratios ):
        rats = ratios.ix[ self.rows, self.cols ]
        ##if np.ndim( rats ) < 2 or np.size( rats, 0 ) <= 1 or np.size( rats, 1 ) <= 1 or \
        ##        np.mean( rats.isnull().values ) > 0.95:
        ##    warnings.warn( "COULD NOT COMPUTE RESIDUE FOR BICLUSTER " + self.k )
        ##    return 1.0
        return funcs.matrix_residue( rats )

    def var( self, ratios, var_add=0.1 ):
        rats = ratios.ix[self.rows, self.cols].copy() ## get the bicluster's submatrix of the data
        mn = rats.mean(0) ## subtract bicluster mean profile
        rats -= mn
        return np.nanvar(rats.values) / (np.nanvar(mn.values) + var_add)

    def re_seed_if_necessary( self, clusters, ratios, min_rows=3, max_rows=80 ):
        ## Need to fill in "empty" clusters, mostly because var is non-defined, Also because meme-ing only works 
        ## with >2 sequences. 
        ## Squash clusters that get way too big (DONE: just remove some to bring it down to max_rows)
        if length(self.rows) < min_rows or length(self.rows) > max_rows:
            nr = length(self.rows)
            warnings.warn( "RESEEDING BICLUSTER %d (%d)" % (self.k, nr) )

        ## DONE: add rows that are preferentially in few (or no) other clusters -- 
        ## sample from get_cluster_row_counts(clusters)
        ##clust.rows = unique([clust.rows, [Distributions.sample([1:nrow(ratios.x)]) for i=1:5]]) ## add 5 random rows
        counts_g = get_all_cluster_row_counts( clusters )
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

    def network_density( self, network ):
        return funcs.subnetwork_density( self.rows, network )

    def meme_pval( self ):
        if np.size(self.mast_out,0) <= 0:
            return NA

        genes = self.rows
        df = self.mast_out[ np.in1d( self.mast_out.Gene, genes ) ]
        mn = np.nanmean( log10( df['P-value'] ) )
        return mn

#     def fill_all_scores( self, iter, ratios, string_net, force=False, verbose=False ):
#         if verbose:
#             print self.k + " " + len(self.rows) + " " + len(self.cols)
#             weight_r, weight_n, weight_m, weight_c, weight_v = get_score_weights(iter)

#             if force or ( ( weight_r + weight_c > 0 ) and ( np.sum(self.changed) > 0 ) ):
#                 try:   ## Actually need to do it for rows too, even if only cols changed
#                     clust = get_cluster_expr_rowcol_scores( clust, ratios )
# ############### CONVERTED UP TO HERE!!!
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


# function copy_cluster( cc::bicluster, full=true, everything=true )
#     out = cc
#     if everything ## make copies of everything
#         out = bicluster( cc.k, copy(cc.rows), copy(cc.cols), cc.var, cc.resid, cc.dens_string, cc.meanp_meme,
#                         copy(cc.scores_r), copy(cc.scores_c), copy(cc.scores_n), copy(cc.scores_m), 
#                         copy(cc.meme_out), copy(cc.mast_out), cc.changed )
#     elseif full ## Make copies of rows/cols so if they get updated in the original, they won't change here.
#         out = bicluster( cc.k, copy(cc.rows), copy(cc.cols), cc.var, cc.resid, cc.dens_string, cc.meanp_meme, 
#                         cc.scores_r, cc.scores_c, cc.scores_n, cc.scores_m, 
#                         cc.meme_out, cc.mast_out, cc.changed )
#     else ## This is a shallow copy - changes to rows/cols in this cluster will mirror those in cc
#         out = bicluster( cc.k, cc.rows, cc.cols, cc.var, cc.resid, cc.dens_string, cc.meanp_meme, 
#                         cc.scores_r, cc.scores_c, cc.scores_n, cc.scores_m, 
#                         cc.meme_out, cc.mast_out, cc.changed )
                   
#     end
#     out
# end
