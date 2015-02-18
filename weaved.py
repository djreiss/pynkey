import scipy.weave

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
