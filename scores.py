## Trial: 1.0*r - 0.9*n + 5.0*m + 0.3*v ... NOTE n has NOT been corrected so we need to negate it here!
## DONE: add volume v0 factor (up-weight to grow clusters)
## DONE: add gene popularity g0 factor (down-weight for adding genes that are already in lots of other clusters and
##      up-weight for genes that are in very few)        
## DONE: do we need weight_c? Might be nice to upweight gains for columns vs. rows or vice versa.
## TODO: Load and use other networks (with their respective scores)

import numpy as np
import pandas as pd

print 'importing scores'

import params
import utils as ut

## Note can make these weights a global var and update them in this function.
## Just need to define them as global IN the function
def get_score_weights(iter, ratios):
    mn = params.max_network_weight
    mm = params.max_motif_weight
    mv = params.max_volume_weight
    mg = params.max_clusters_per_gene_weight
    mc = params.max_column_weight
    nit = params.n_iters
    weight_r =   1.0
    weight_n =   0.0 + -mn * float(iter-1) / nit   ## increase linearly from 0 at iter=1 to 0.9
    weight_m =  (1.0 +  mm * float(iter-1) / nit) * (0 if iter<=5 else 1) ## ramp up from 1 to 1.8 starting at iter=6
    ## GOOD: weight_n =   0.5 + -mn * float(iter-1) / nit   ## increase linearly from 0 at iter=1 to 0.9
    ## GOOD: weight_m =  (0.5 +  mm * float(iter-1) / nit) ##* (0 if iter<=5 else 1) ## ramp up from 1 to 1.8 starting at iter=6
    weight_c =   0.0 +  mc * np.size(ratios,0)/np.size(ratios,1)/1.2 ## ??? ## 1.2 works good for hpy
    weight_v =   0.1 +  mv * float(iter-1) / nit  ## ramp up from 0.3 to 0.8
    weight_g =   0.1 +  mg * float(iter-1) / nit  ## ramp up from 0.3 to 0.8
    return (weight_r, weight_n, weight_m, weight_c, weight_v, weight_g)

## TODO: use numexpr to speed up and avoid temporary array creation
## also use array.fill(0) to reset to zero without recreating a new array
## scores_DF has columns ['score_r', 'score_n', 'score_m', 'score_c', 'score_v', 'score_g']
def get_combined_scores( scores_DF, iter, ratios, do_stdize=True ): 
    weight_r, weight_n, weight_m, weight_c, weight_v, weight_g = get_score_weights( iter, ratios )
    wts = pd.Series( [weight_r, weight_n, weight_m, weight_c, weight_v, weight_g], 
                     index=['score_r', 'score_n', 'score_m', 'score_c', 'score_v', 'score_g'] )

    df = scores_DF
    nr = df.shape[0]
    out = np.zeros( nr )
    tmp = np.zeros_like( scores_DF.score_n )

    for name in ['score_r', 'score_m', 'score_n', 'score_v', 'score_g']:  ## TBD: what about score_c ??
        wt = wts[name]
        if abs(wt) < 1e-5:
            continue
        scor = df[name].values
        isnans = np.isnan(scor)
        if np.all(isnans):
            continue
        if do_stdize:     ## Need to standardize each score first
            scor = ut.sdize_vector(scor)
        out[ ~isnans ] += wt * scor[ ~isnans ]

    if do_stdize:
        out = ut.sdize_vector(out)
    return out

def get_combined_score( r, n, m, v, g, iter, ratios ):
    weight_r, weight_n, weight_m, weight_c, weight_v, weight_g = get_score_weights( iter, ratios )

    def zero_if_bad(x):
        if np.isnan(x) or np.isinf(x):
            return 0
        return x

    out = weight_r * zero_if_bad(r) + \
        weight_n * zero_if_bad(n) + \
        weight_m * zero_if_bad(m) + \
        weight_v * zero_if_bad(v) + \
        weight_g * zero_if_bad(g)
    return out

