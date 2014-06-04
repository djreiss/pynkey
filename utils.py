from numpy.random import uniform
import numpy as np
import random
 
def slice_sampler(px, N = 1, x = None):
    """
    Provides samples from a user-defined distribution.
    
    slice_sampler(px, N = 1, x = None)
    
    Inputs:
    px = A discrete probability distribution.
    N  = Number of samples to return, default is 1
    x  = Optional list/array of observation values to return, where prob(x) = px.
 
    Outputs:
    If x=None (default) or if len(x) != len(px), it will return an array of integers
    between 0 and len(px)-1. If x is supplied, it will return the
    samples from x according to the distribution px.    

    From:
    http://www.adamlaiacano.com/post/14987215771/python-function-for-sampling-from-an-arbitrary-discrete
    """
    values = np.zeros(N, dtype=numpy.int)
    samples = np.arange(len(px))
    px = np.array(px) / (1.*sum(px))
    u = uniform(0, max(px))
    for n in xrange(N):
        included = px>=u
        choice = random.sample(range(np.sum(included)), 1)[0]
        values[n] = samples[included][choice]
        u = uniform(0, px[included][choice])
    if x:
        if len(x) == len(px):
            x=np.array(x)
            values = x[values]
        else:
            print "px and x are different lengths. Returning index locations for px."
    if N == 1:
        return values[0]
    return values

def sdize_vector( vec ):
    v = vec[ vec != 0 ]
    mn = np.nanmean( v )
    sd = np.nanstd( v )
    return (vec - mn) / ( sd + 0.001 )

def minmax( vec ):
    return (np.nanmin(vec), np.nanmax(vec))
