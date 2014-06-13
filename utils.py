import tempfile

from numpy.random import uniform
import numpy as np
import random

import multiprocessing as mp
## see: http://matthewrocklin.com/blog/work/2013/12/05/Parallelism-and-Serialization/
## and https://stackoverflow.com/questions/1816958/cant-pickle-type-instancemethod-when-using-pythons-multiprocessing-pool-ma
## and https://stackoverflow.com/questions/19984152/what-can-multiprocessing-and-dill-do-together
##import pathos.multiprocessing as mp

## Note these par funcs need to be in globals.py (and hence global) to have global access to all
##    data -- I don't know how to change this right now.

print 'importing utils'

def do_something_par( items, func, threads=None ): # None ensures that it makes as many threads as avail. in computer
    if threads == 1:
        out = map(func, items)
    else:
        pool = mp.Pool(processes=threads)              # start 4 worker processes
        ## Need to send, e.g. a tuple (1, counts_g) if fill_all_scores_par() took multiple args
        out = pool.map(func, items)
        pool.terminate()
    return out

def writeLines( lines, fname=None ):
    if fname is None:
        fname = tempfile.mktemp()
    handle = open(fname, 'w')
    handle.write('\n'.join(lines)) ## assumes lines are an array split by '\n' - if not then do '\n'.join(lines) first
    handle.close()
    return fname
 
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
    values = np.zeros(N, dtype=np.int)
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

def sdize_vector( vec, ignore_zeroes=True ): ## note this is inplace! If don't want, pass vec.copy() !!
    v = vec
    if ignore_zeroes:
        v = vec[ vec != 0 ]
    mn = np.nanmean( v )
    sd = np.nanstd( v )
    vec -= mn
    vec /= (sd + 0.001) ## try to minimize copies?
    return vec

def minmax( vec ):
    return (np.nanmin(vec), np.nanmax(vec))
