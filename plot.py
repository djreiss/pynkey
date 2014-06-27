import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

print 'importing plot'

import globals
import params

def setup_text_plots(fontsize=8, usetex=True):
    """
    This function adjusts matplotlib settings so that all figures in the
    textbook have a uniform format and look.
    """
    import matplotlib
    matplotlib.rc('legend', fontsize=fontsize, handlelength=3)
    matplotlib.rc('axes', titlesize=fontsize, labelsize=fontsize)
    matplotlib.rc('xtick', labelsize=fontsize)
    matplotlib.rc('ytick', labelsize=fontsize)
    matplotlib.rc('text', usetex=usetex)
    ## see https://stackoverflow.com/questions/12322738/how-do-i-change-the-axis-tick-font-in-a-matplotlib-plot-when-rendering-using-lat
    ##matplotlib.rc('text.latex', preamble=r'\usepackage{cmbright}')
    matplotlib.rc('mathtext', fontset='stixsans')
    matplotlib.rc('font', size=fontsize, family='sans-serif',
                  style='normal', variant='normal',
                  stretch='normal', weight='normal')

def plot_stats():
    stats = globals.stats_df
    print stats.tail()

## This is a good intro to pandas dataframes
## http://pandas.pydata.org/pandas-docs/stable/10min.html

## try plotting it. using astroML/book_figures/chapter1/fig_SDSS_imaging.py as an example
## this fuinction sets up the plots to look the same as in the astroML text
        
    setup_text_plots(fontsize=8, usetex=True)

    ##plot_kwargs = dict(color='k', linestyle='none', marker=',')
    plot_kwargs = dict(color='k', marker=',')
    plt.close() ## close previous window if exists
    fig = plt.figure(figsize=(10, 8))

    # Now plot with matplotlib 
    ax1 = fig.add_subplot(3,3,1)
    ax1.plot(stats.iter, stats.RESID, **plot_kwargs)
    ax1.set_ylabel('Residual')
    ax1.set_xlabel('Iteration')

    ax2 = fig.add_subplot(3,3,2) ##, sharex=ax1)
    vals = stats[["iter", "MEME_PVAL"]].dropna(how="any").values
    ax2.plot(vals[:, 0], vals[:, 1], **plot_kwargs)
    ax2.set_ylabel('Motif p-value')
    ax2.set_xlabel('Iteration')

    ax3 = fig.add_subplot(3,3,3)
    vals = stats[["iter", "STRING_DENS"]].dropna(how="any").values
    ax3.plot(vals[:, 0], vals[:, 1], **plot_kwargs)
    ax3.set_ylabel('Avg STRING Network Density')
    ax3.set_xlabel('Iteration')

    ax4 = fig.add_subplot(3,3,4)
    plot_kwargs = dict(marker=',')
    vals = stats[["iter", "r0"]].dropna(how="any").values
    ax4.plot(vals[:, 0], vals[:, 1], color='k', **plot_kwargs)
    vals = stats[["iter", "m0"]].dropna(how="any").values
    ax4.plot(vals[:, 0], vals[:, 1], color='r', **plot_kwargs)
    vals = stats[["iter", "n0"]].dropna(how="any").values
    ax4.plot(vals[:, 0], vals[:, 1], color='g', **plot_kwargs)
    vals = stats[["iter", "c0"]].dropna(how="any").values
    ax4.plot(vals[:, 0], vals[:, 1]/10.0, color='b', **plot_kwargs)
    vals = stats[["iter", "v0"]].dropna(how="any").values
    ax4.plot(vals[:, 0], vals[:, 1], color='c', **plot_kwargs)
    vals = stats[["iter", "g0"]].dropna(how="any").values
    ax4.plot(vals[:, 0], vals[:, 1], color='m', **plot_kwargs)
    ax4.set_xlabel('Iteration')
    ax4.set_ylabel('Scaling')

    text_kwargs = dict(transform=plt.gca().transAxes, ha='center', va='bottom')
    ax4.text(0.12, 0.9, 'r0', color='k', **text_kwargs)
    ax4.text(0.24, 0.9, 'm0', color='r', **text_kwargs)
    ax4.text(0.36, 0.9, 'n0', color='g', **text_kwargs)
    ax4.text(0.48, 0.9, 'c0/10', color='b', **text_kwargs)
    ax4.text(0.60, 0.9, 'v0', color='c', **text_kwargs)
    ax4.text(0.72, 0.9, 'g0', color='m', **text_kwargs)
    
    ax2a = fig.add_subplot(3,3,5) ##, sharex=ax1)
    vals = stats[["iter", "ROWS"]].dropna(how="any").values
    ax2a.plot(vals[:, 0], vals[:, 1], color='b', **plot_kwargs)
    vals = stats[["iter", "COLS"]].dropna(how="any").values
    ax2a.plot(vals[:, 0], vals[:, 1], color='r', **plot_kwargs)
    ax2a.set_ylabel('Avg Number')
    text_kwargs = dict(transform=plt.gca().transAxes, ha='center', va='bottom')
    ax2a.text(0.25, 0.9, 'Rows', color='b', **text_kwargs)
    ax2a.text(0.75, 0.9, 'Cols', color='r', **text_kwargs)
    ax2a.set_xlabel('Iteration')
    
    ## Plot histograms of #genes, #conds in each bicluster
    ax4a = fig.add_subplot(3,3,6)
    nrows = np.array( [len(c.rows) for c in globals.clusters.values()] )
    plt.hist( nrows, 20 )
    ax4a.set_xlabel('Number per cluster')
    ax4a.set_ylabel('Count')

    ncols = np.array( [len(c.cols) for c in globals.clusters.values()] )
    plt.hist( ncols, 20, color='r' )
    text_kwargs = dict(transform=plt.gca().transAxes,
                       ha='center', va='bottom')
    ax4a.text(0.20, 0.9, 'Rows', color='b', **text_kwargs)
    ax4a.text(0.80, 0.9, 'Cols', color='r', **text_kwargs)

    ax5 = fig.add_subplot(3,3,7)
    vals = stats[["iter", "N_MOVES"]].dropna(how="any").values
    ax5.plot(vals[:, 0], vals[:, 1], color='k', **plot_kwargs)
    vals = stats[["iter", "N_IMPROVEMENTS"]].dropna(how="any").values
    ax5.plot(vals[:, 0], vals[:, 1], color='r', **plot_kwargs)
    text_kwargs = dict(transform=plt.gca().transAxes,
                       ha='center', va='bottom')
    ax5.text(0.20, 0.9, 'Num Moves', color='k', **text_kwargs)
    ax5.text(0.70, 0.9, 'Num Improvements', color='r', **text_kwargs)
    ax5.set_xlabel('Iteration')
    ax5.set_ylabel('Count')

    ax5a = fig.add_subplot(3,3,8)
    resids = np.array( [c.resid for c in globals.clusters.values()] )
    resids = resids[ resids > 0.001 ]
    if len(resids) > 0:
        plt.hist( resids, 20 )
        ax5a.set_xlabel('Cluster Residuals')
        ax5a.set_ylabel('Count')

    ax6 = fig.add_subplot(3,3,9)
    pclusts = np.array( [c.meanp_meme for c in globals.clusters.values()] )
    pclusts = pclusts[ ~np.isnan(pclusts) & ~np.isinf(pclusts) ]
    if len(pclusts) > 0:
        plt.hist( pclusts, 20 )
        ax6.set_xlabel('Cluster Mean log10(P-value)s')
        ax6.set_ylabel('Count')

    plt.show()

# def plot_scores(scores): ## plot scores dataframe from floc.get_scores_all()
#     setup_text_plots(fontsize=8, usetex=True)

#     ##plot_kwargs = dict(color='k', linestyle='none', marker=',')
#     plot_kwargs = dict(color='k', marker=',')
#     plt.close() ## close previous window if exists
#     fig = plt.figure(figsize=(10, 8))

#     ax1 = fig.add_subplot(3,3,1)
#     ax1.plot(stats.iter, stats.RESID, **plot_kwargs)
#     ax1.set_ylabel('Residual')
#     ax1.set_xlabel('Iteration')

def plot_motif_logo( memeOut, motInd=1 ):
    import meme
    import matplotlib.image as mpimg

    record = meme.parseMemeOut( memeOut )
    kwargs = dict(color_scheme='classic')
    record[motInd].weblogo('file.png', color_scheme='color_classic') ## note, can use format='PDF'
    img = mpimg.imread('file.png')
    imgplot = plt.imshow( img )

if __name__ == '__main__':
    ## see https://docs.python.org/2/library/optparse.html
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-o", "--organism", dest="organism",
                  help="input organism [default: %default]", metavar="ORGANISM", default='Hpy')
    #parser.add_option("-i", "--iter", dest="iter", type=int,
    #              help="desired iteration to plot [default: %default]", metavar="ITER", default='100')
    (options, args) = parser.parse_args()

    import init
    import os.path
    for iter in xrange(1,100):
        if os.path.exists( 'output/%s_%04d.pkl' % (options.organism, iter+1) ):
            continue
        else:
            init.init( 'output/%s_%04d.pkl' % (options.organism, iter) )
            break
    
    plot_stats()
