## Trial: 1.0*r - 0.9*n + 5.0*m + 0.3*v ... NOTE n has NOT been corrected so we need to negate it here!
## DONE: add volume v0 factor (up-weight to grow clusters)
## DONE: add gene popularity g0 factor (down-weight for adding genes that are already in lots of other clusters and
##      up-weight for genes that are in very few)        
## DONE: do we need weight_c? Might be nice to upweight gains for columns vs. rows or vice versa.
## TODO: Load and use other networks (with their respective scores)

import numpy as np

import params

## Note can make these weights a global var and update them in this function.
## Just need to define them as global IN the function
def get_score_weights(iter, ratios):
    ##global n_iters, ratios, max_network_weight, max_motif_weight
    mn = params.max_network_weight
    mm = params.max_motif_weight
    nit = params.n_iters
    weight_r =  1.0
    weight_n =  0.0 + -mn * float(iter-1) / nit   ## increase linearly from 0 at iter=1 to 0.9
    weight_m =  1.0 +  mm * float(iter-1) / nit * (0 if iter<=5 else 1) ## ramp up from 1 to 1.8 starting at iter=6
    weight_c =  0.0 + 0.2 * np.size(ratios,0)/np.size(ratios,1)/12.0 ## ??? ## 1.2 works good for hpy
    weight_v =  0.1 + 0.3 * float(iter-1) / nit  ## ramp up from 0.3 to 0.8
    weight_g =  0.1 + 0.1 * float(iter-1) / nit  ## ramp up from 0.3 to 0.8
    return weight_r, weight_n, weight_m, weight_c, weight_v, weight_g

## TODO: use numexpr to speed up and avoid temporary array creation
## also use array.fill(0) to reset to zero without recreating a new array
## scores_DF has columns ['score_r', 'score_n', 'score_m', 'score_c', 'score_v', 'score_g']
def get_combined_scores( scores_DF, iter, ratios ): 
    weight_r, weight_n, weight_m, weight_c, weight_v, weight_g = get_score_weights( iter, ratios )

    df = scores_DF
    nr = df.shape[0]
    ## Need to standardize each scores column first
    out = np.zeros( nr )
    out[ np.invert(np.isnan(df.score_r.values)) ] += weight_r * df.score_r[ np.invert(np.isnan(df.score_r.values)) ]
    tmp = np.zeros_like( scores_DF.score_n )
    tmp[ np.invert(np.isnan(df.score_n.values)) ] += weight_n * df.score_n[ np.invert(np.isnan(df.score_n.values)) ]
    out += tmp
    tmp[:] = 0
    tmp[ np.invert(np.isnan(df.score_m.values)) ] += weight_m * df.score_m[ np.invert(np.isnan(df.score_m.values)) ]
    out += tmp
    tmp[:] = 0
    tmp[ np.invert(np.isnan(df.score_v.values)) ] += weight_v * df.score_v[ np.invert(np.isnan(df.score_v.values)) ]
    out += tmp
    tmp[:] = 0
    tmp[ np.invert(np.isnan(df.score_g.values)) ] += weight_g * df.score_g[ np.invert(np.isnan(df.score_g.values)) ]
    out += tmp
    
# function get_combined_score( r::Float32, n::Float32, m::Float32, v::Float32, g::Float32 )
#     (weight_r, weight_n, weight_m, weight_c, weight_v, weight_g) = get_score_weights()
#     weight_r * ( isnan(r) || isinf(r) ? 0 : r ) + 
#     weight_n * ( isnan(n) || isinf(n) ? 0 : n ) + 
#     weight_m * ( isnan(m) || isinf(m) ? 0 : m ) +
#     weight_v * ( isnan(v) || isinf(v) ? 0 : v ) +
#     weight_g * ( isnan(g) || isinf(g) ? 0 : g )
# end
