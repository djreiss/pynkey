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

    def compute_residue( self, ratios ):
        rats = ratios.ix[ self.rows, self.cols ]
        ##if np.ndim( rats ) < 2 or np.size( rats, 0 ) <= 1 or np.size( rats, 1 ) <= 1 or \
        ##        np.mean( rats.isnull().values ) > 0.95:
        ##    warnings.warn( "COULD NOT COMPUTE RESIDUE FOR BICLUSTER " + self.k )
        ##    return 1.0
        return funcs.matrix_residue( rats )

    def compute_residue_deltas( self, ratios, all_genes ):
        is_in = np.in1d( all_genes, self.rows )
        rows = self.rows
        rats = ratios.ix[ :, self.cols ]
        resid = funcs.matrix_residue( ratios.ix[ rows ] )
        all_resids = np.zeros( len( all_genes ), float )
        for i in xrange(len(all_genes)):
            r = all_genes[i]
            rows2 = np.append( rows, r ) if is_in[i] else rows[ rows != r ]
            all_resids[i] = funcs.matrix_residue( rats.ix[ rows2, ] )
        return all_resids - resid

    def compute_var( self, ratios, var_add=0.1 ):
        rats = ratios.ix[self.rows, self.cols] ##.copy() ## get the bicluster's submatrix of the data
        ##mn = nanmean(rats,0) ## subtract bicluster mean profile
        ##rats -= mn
        ##return np.nanvar(rats.values) / (np.nanvar(mn.values) + var_add)
        return funcs.matrix_var( rats, var_add )

    def compute_network_density( self, network ):
        return funcs.subnetwork_density( self.rows, network )

    def compute_network_density_deltas( self, network, all_genes ):
        is_in = np.in1d( all_genes, self.rows )
        rows = self.rows
        net1 = network.ix[ self.rows ]
        dens = funcs.subnetwork_density( rows, net1 )
        all_dens = np.zeros( len( all_genes ), float )
        for i in xrange(len(all_genes)):
            r = all_genes[i]
            rows2 = np.append( rows, r ) if is_in[i] else rows[ rows != r ]
            all_dens[i] = funcs.subnetwork_density( rows2, net1 )
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
    def get_volume_row_scores( self, all_genes ): ##, is_in_r )
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
    def get_volume_col_scores( self, ratios ): ##, is_in_c )
        is_in = np.in1d( ratios.cols.values, self.cols )
        lc = len(self.cols)
        score_vc = np.array( [ +1.0/lc if i else -1.0/lc for i in is_in ] )
        ##score_vc = fill(NA, length(b.scores_c)) ## Do we really care how many cols there are?
        return score_vc

    ## Up-weight moves OUT if counts_g is HIGH, and moves IN if counts_g is LOW
    ## counts_g comes from counts_g = funcs.get_all_cluster_row_counts( clusters, all_genes )
    def get_row_count_scores( self, counts_g, all_genes ):
        is_in = np.in1d( all_genes, self.rows )
        score_g = get_cluster_row_count_scores( counts_g )
        score_g = np.array( [ (+score_g[all_genes[i]] if is_in[i] else -score_g[all_genes[i]]) for i in xrange(len(is_in)) ] )
        return score_g

    ## counts_g comes from counts_g = funcs.get_all_cluster_row_counts( clusters, all_genes )
    def fill_all_scores(self, counts_g, all_genes, ratios, string_net):
        self.scores_r = self.compute_residue_deltas( ratios, all_genes )
        self.scores_n = self.compute_network_density_deltas( string_net, all_genes )

    ## counts_g comes from counts_g = funcs.get_all_cluster_row_counts( clusters, all_genes )
    def get_floc_scoresDF_rows(self, counts_g, all_genes, iter, ratios):
        is_in = np.in1d( all_genes, self.rows)
        ## only use them if their weights are > 0
        weight_r, weight_n, weight_m, weight_c, weight_v, weight_g = scores.get_score_weights( iter, ratios )
        NAs =  np.repeat(NA, len(self.scores_r)) if (weight_r <= 0 or weight_n <= 0 or weight_m <= 0) else NA
        score_r = self.scores_r if weight_r > 0 else NAs
        score_n = self.scores_n if abs(weight_n) > 0 else NAs
        score_m = self.scores_m if weight_m > 0 else NAs
        score_vr = self.get_volume_row_scores( all_genes )
        score_g = self.get_row_count_scores( counts_g, all_genes )
        out = pd.DataFrame( { 'row_col_ind':all_genes,
                              'is_in':is_in,
                              'is_row_col':np.repeat('r', len(self.scores_r)), ## CANT: move this outside the loop
                              'k':np.repeat(self.k, len(self.scores_r)),
                              'score':score_r,
                              'score_n':score_n,
                              'score_m':score_m,
                              'score_v':score_vr,
                              'score_g':score_g } )
        return out

#     def get_expr_rowcol_scores( self, ratios, var_add=0.1 ):
#      ## Try getting the effect of adding/removing each row/col from this cluster on its total variance
#      ## DONE: want to have biclusters with larger row-variance, so normalize var(x)/var(mn)
#      ## Other variables to be tweaked: currently using out_var(x)=var(x)/(var(colmeans(x))+1) --
#      ##       perhaps the "1" should be lower... see "var_add" variable
#      ## Commented out volume factor - shoul add this in to total, but externally to this function
#         rows=self.rows; cols=self.cols; score_r=self.scores_r; score_c=self.scores_c

#      ##const bsize = length(rows)*length(cols)
#         vr = self.compute_var( ratios )
#         resid = self.compute_residue( ratios )

#         rats = ratios.ix[self.rows, self.cols].copy() ## get the bicluster's submatrix of the data
#         mn = nanmean(rats,0) ## subtract bicluster mean profile
#         rats -= mn

#      ##xx::Matrix{T} = x.x[:,cols] ## This uses less RAM than previous version; about same speed
#         xx = ratios.ix[:,cols]
#         mn = np.nanmean(xx.ix[rows,:], 0)
#         x2 = xx.copy()  ## Get the full matrix; subtract the mean biclust profile
#         x2.apply( lambda x: x-mn )
    
#      ## Get changes in variance for adding/removing rows; normalize by change in cluster volume. Lower is better!
#      ##v_factor::T = length(cols) / bsize
#         if len(self.scores_r) == 0:
#             self.scores_r = np.empty(size(ratios,0))
# #     tmpVec::Vector{T} = Array(T, size(xx,2))
#         for r in ratios.index.values: ## Iterate over rows
#             isIn = r in rows ## r is in the bicluster
#             newR = np.delete(rows,np.where(rows==r)) if isIn else np.append(rows,r) ## TBD: prevent reallocation of vector here?
#             newvr = np.nanvar(x2.ix[newR,:]) / ( np.nanvar(np.nanmean(xx.ix[newR,:],1) ) ) + var_add ) ## DONE -- prevent realloc of matrix here
#        ## NOTE: the colmeans() above is the slow part here -- convince yourself that you need it!
#        score_r[r] = ( newvr - vr ) #/ (vr+0.01) ##+ ( isIn ? v_factor : -v_factor )
#     end

#     ##xx=x.x[rows,:]
#     xx=view(x.x,rows,[1:size(x.x,2)])
#     mn=colmeans(xx)
#     x2=similar(xx) ## Get the full matrix; subtract the mean biclust profile
#     for i=1:size(x2,2) x2[:,i] = xx[:,i] - mn[i]; end

#     ##v_factor = length(rows) / bsize
#     if length(score_c) == 0 score_c = Array(T,size(x.x,2)); end
#     for c in 1:size(x.x,2) ## Iterate over cols
#        isIn = in(cols, c) ##any(cols .== c) ## c is in the bicluster
#        newC = isIn ? remove(cols,c) : [cols, [c]] ##append(cols, c) ##[cols,c]
#        newvr = nanvar(x2[:,newC]) / ( nanvar(mn[newC]) + var_add )
#        score_c[c] = ( newvr - vr ) #/ (vr+0.01) ##+ ( isIn ? v_factor : -v_factor )
#     end

#     b.var = vr
#     b.resid = resid
#     b.scores_r = score_r ## Not needed since the vector is filled (passed by reference)
#     b.scores_c = score_c ## but we'll do it here anyway b/c it doesn't slow things down much
#     b
# end

# #@iprofile begin
# ## NOTE WE ASSUME THE NETWORK IS "SYMMETRIZED" if it is undirected - i.e. (n1,n2,w) and (n2,n1,2) are both in the dataframe.
# ## NOTE we return log10(density) changes, which is OPPOSITE of ROW/MOTIF scores because increases are better
# function get_cluster_network_row_scores( b::bicluster, network::DataFrame ) 
#     global ratios, all_genes
#     ## Note that right now rows in the bicluster is the index into the rownames of the expression matrix!!! TODO: change this!
#     r_rownames = rownames(ratios)
#     rows = r_rownames[b.rows] ##keys( all_genes )[ in( values( all_genes ), b.rows ) ]
#     net = sub( network, findin( network["x1"], rows ) ) ## Assume network is symmetric!
#     grps = groupby( net, "x2" ) ## each group is a subnetwork with all genes in bicluster (x1) connected to a given gene
#     grpnames = convert( Vector{ASCIIString}, [ grps[i]["x2"][1] for i=1:length( grps ) ] ) ## the given gene for each group
#     grpname_lookup = [ grpnames[i] => i for i=1:length(grps) ]
#     new_net = net2 = sub( net, findin( net["x2"], rows ) )   ## subnetwork just connecting between nodes in the bicluster
#     sum_weights = sum( new_net["x3"] )
#     dens = log10( sum_weights / length(rows)^2 + 1e-9 ) ## Already symmetrized, need to decrease count by 1/2
#     n2 = net2["x2"] ##Vector{ASCIIString} = net2[2]

#     score_n = b.scores_n
#     if length(score_n) == 0 score_n = Array(Float32,size(ratios.x,1)); end    
#     new_dens = 0.
#     for r in r_rownames
#         isIn = in(rows, r) 
#         if isIn ## r is in the bicluster -- remove the node from the bicluster-only subnetwork and recalc the density
#             newR = remove(rows,r)
#             new_net = sub( net2, findin( n2, newR ) ) ##net2[2].data, newR ) ) )
#             new_dens = sum( new_net["x3"] ) / length(newR) / length(rows)
#         else
#             try 
#                 new_net = grps[ grpname_lookup[r] ]
#                 new_dens = ( sum_weights + sum( new_net["x3"] ) ) / (length(rows)+1) / length(rows)
#             catch ## no subnetwork for this "r"; so use the bicluster's subnetwork, but diminished by increased length(newR)
#                 new_dens = sum_weights / (length(rows)+1) / length(rows)
#             end
#         end
#         score_n[all_genes[r]] = log10(new_dens+1e-9) - dens ##/ (dens+0.1) ##+ ( isIn ? v_factor : -v_factor )  ## lower is BETTER
#     end
#     b.scores_n = score_n
#     b.dens_string = dens
#     b
# end
# #end # profiler

# function get_cluster_meme_row_scores( b::bicluster )
#     global ratios, all_genes
#     #println("IN HERE: MOT SCORES $(b.k)")

#     score_m = b.scores_m
#     if length(score_m) == 0 score_m = Array(Float32,size(ratios.x,1)); end
#     r_rownames = rownames(ratios)
#     genes = r_rownames[b.rows]
#     #sz=size(b.mast_out);#println("HERE MOT1 $sz")
#     ## TODO: don't only use values for genes that are in the ratios matrix.
#     ## DONE: make mast_out into a DataFrame for easier searching, subsetting
#     if size(b.mast_out,1) <= 0
#         b.scores_m = float32( zeros( size(ratios.x,1) ) )
#         b.meanp_meme = NA
#         return( b )
#     end
#     df = sub( b.mast_out, findin( b.mast_out["Gene"], genes ) )
#     #sz=size(df);#println("HERE MOT2 $sz")
#     mn::Float32 = nanmean( log10( df["P-value"].data ) )
#     pvals::Dict{ASCIIString,Float32} = Dict{ASCIIString,Float32}()
#     for i in 1:size(b.mast_out,1) pvals[b.mast_out["Gene"].data[i]] = b.mast_out["P-value"].data[i]; end
#     not_there = r_rownames[ ! in(r_rownames, b.mast_out["Gene"].data) ]
#     for r in not_there pvals[r] = NA; end
#     for g in r_rownames ## Iterate over rows
#        isIn = in(genes, g) ##any(rows .== r) ## r is in the bicluster
#        newR = isIn ? remove(genes,g) : [genes,[g]]
#        pvs = float32( [ pvals[r] for r in newR ] )
#        # pvs = df["P-value"].data
#        # if isIn pvs = sub( df, df["Gene"].data .!= g )["P-value"].data
#        # else    pvs = [ pvs, sub( b.mast_out, b.mast_out["Gene"].data .== g )["P-value"].data ]
#        # end
#        newmn = nanmean( log10( pvs ) )
#        score_m[all_genes[g]] = newmn - mn ##/ (mn+0.1) ##+ ( isIn ? v_factor : -v_factor )
#     end

#     b.scores_m = score_m
#     b.meanp_meme = mn
#     b
# end

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
        counts_g = funcs.get_all_cluster_row_counts( clusters, all_genes )
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


## THis is a static function so outside the definition of the bicluster
## counts_g comes from counts_g = funcs.get_all_cluster_row_counts( clusters, all_genes )
def get_cluster_row_count_scores( counts_g ):
    thresh = params.avg_clusters_per_gene ##1.3 ## 2.0 ## 3.0 ## lower is better; coerce removing if gene is in more than 2 clusters
    return dict( { (i,j-thresh) for i,j in counts_g.items() } )

## TBD: write a generic update scores function that tests "update_func" (e.g. compute_resid) if a
##    each gene is added/removed from the cluster
def get_update_scores( cluster, update_func, *args ):
    print args
    return update_func( cluster, *args )
