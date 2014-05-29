import sys

from numpy import nan as NA
import numpy as np
from numpy import random as rand
import pandas as pd

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

        max_float = sys.float_info.max
        self.var = None ##max_float
        self.resid = None ##max_float
        self.dens_string = None ##max_float
        self.meanp_meme = None ##max_float
        self.scores_r = np.array([])
        self.scores_c = np.array([])
        self.scores_n = np.array([])
        self.scores_m = np.array([])
        self.meme_out = ['']
        self.mast_out = pd.DataFrame()
        self.changed = np.ones(2, np.bool)

    def __repr__(self):
        return 'Bicluster: %d' % self.k + '\n' + \
            'resid: %f' % self.resid + '\n' + \
            'meme-p: %f' % self.meanp_meme + '\n' + \
            'string: %f' % self.dens_string + '\n' + \
            'rows: %s' % str(self.rows) + '\n' + \
            'cols: %s' % str(self.cols)
    
