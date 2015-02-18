import numpy as np
import scipy.weave
from numpy import random as rnd

import utils as ut

print 'importing weaved'

def fast_resid( rats ):
    code = """
    blitz::Array<double,1> d_rows(Nrats[0]), d_cols(Nrats[1]);
    blitz::Array<double,1> n_cols(Nrats[1]);
    d_cols = 0.0; n_cols = 0.0; 
    double d_all = 0, n_all = 0;
    int nr = Nrats[0], nc = Nrats[1];
    for ( int i=0; i<nr; i++ ) {
      d_rows(i) = 0.0; // can initialize these in the loop
      double n_rows_i = 0.0; 
      for ( int j=0; j<nc; j++ ) {
        double rr = rats(i,j);
        if ( isnan(rr) ) continue;
        d_rows(i) += rr;
        d_cols(j) += rr;
        d_all += rr;
        n_rows_i ++;
        n_cols(j) ++;
        n_all ++;
      }
      d_rows(i) /= n_rows_i;
    }
    d_all /= n_all;
    double resid = 0;
    for ( int j=0; j<nc; j++ ) {
      d_cols(j) /= n_cols(j);
      for ( int i=0; i<nr; i++ ) {
        double rr = rats(i,j);
        if ( isnan(rr) ) continue;
        resid += fabs(rr - d_rows(i) - d_cols(j) + d_all); //, 2.0);
      }
    }
    resid /= n_all;
    return_val = resid;
"""
    #d_rows = np.zeros( rats.shape[0] )
    #d_cols = np.zeros( rats.shape[1] )
    #n_rows = np.zeros( rats.shape[0] )
    #n_cols = np.zeros( rats.shape[1] )
    ## See here for other weave options:
    ## https://mail.python.org/pipermail/cplusplus-sig/2002-January/000428.html
    resid = scipy.weave.inline(code, ['rats'], ##, 'd_rows', 'd_cols'], ##, 'n_rows', 'n_cols'],
                   type_converters=scipy.weave.converters.blitz,
                   compiler='gcc', extra_compile_args=['-O3','-malign-double','-funroll-loops'],
                   headers=['<math.h>', '<blitz/array.h>'])
    return resid



## this is weaved (using c++) -- takes about 3.04 seconds
def rnd_bubblesort3( scores, Nrepeats=None ): ## make sure scores is a copy, b/c NaNs will get replaced in this copy
    lsc = len(scores)
    if Nrepeats is None:
        Nrepeats = lsc * 2
    ords = np.arange(lsc)
    rnd.shuffle(ords) ## start w/ random order
    tmp = ut.minmax(scores)
    R = 2.0 * ( tmp[1]-tmp[0] ) ## Denominator of value to compute prob. from
    the_max = tmp[1]
    n = lsc - 1
    sc = scores.copy()
    sc[ np.isnan(sc) ] = the_max ## replace NaN with maximum score
    switchesN = np.array([0]) ## count the number of switches. Not really necessary
    for i in xrange(Nrepeats):
        rnds = rnd.rand(n)
        code = """
        for ( int j=0; j<n; j++ ) {
          int o1 = ords(j), o2 = ords(j+1);
          double g1 = sc(o1), g2 = sc(o2);
          if ( g1 == g2 && g2 == (double)the_max ) continue;
          double p = 0.5 + ( g1 - g2 ) / (double)R; // compute prob of switching
          if ( rnds(j) < p ) { // switch???
            ords(j) = o2;
            ords(j+1) = o1;
            switchesN(0) ++;
          }
        }
        """
        ## See here for other weave options:
        ## https://mail.python.org/pipermail/cplusplus-sig/2002-January/000428.html
        scipy.weave.inline(code,['ords','n','sc','R','switchesN','rnds','the_max'],
                           type_converters=scipy.weave.converters.blitz,
                           compiler='gcc', extra_compile_args=['-O3','-malign-double','-funroll-loops'])
        if i % 1000 == 1:
            print i, switchesN[0], Nrepeats
    return ords
