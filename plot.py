import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

import globals

def plot_stats():
    stats = globals.stats_df
    stats.head()

## This is a good intro to pandas dataframes
## http://pandas.pydata.org/pandas-docs/stable/10min.html

## try plotting it. using astroML/book_figures/chapter1/fig_SDSS_imaging.py as an example
## this fuinction sets up the plots to look the same as in the astroML text

    def setup_text_plots(fontsize=8, usetex=True):
        """
        This function adjusts matplotlib settings so that all figures in the
        textbook have a uniform format and look.
        """
        import matplotlib
        matplotlib.rc('legend', fontsize=fontsize, handlelength=3)
        matplotlib.rc('axes', titlesize=fontsize)
        matplotlib.rc('axes', labelsize=fontsize)
        matplotlib.rc('xtick', labelsize=fontsize)
        matplotlib.rc('ytick', labelsize=fontsize)
        matplotlib.rc('text', usetex=usetex)
        matplotlib.rc('font', size=fontsize, family='serif',
                      style='normal', variant='normal',
                      stretch='normal', weight='normal')
        
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

    ax2 = fig.add_subplot(3,3,2, sharex=ax1)
    vals = stats[["iter", "MEME_PVAL"]].dropna(how="any").values
    ax2.plot(vals[:, 0], vals[:, 1], **plot_kwargs)
    ax2.set_ylabel('Motif p-value')
    ax2.set_xlabel('Iteration')

    ax3 = fig.add_subplot(3,3,3, sharex=ax1)
    vals = stats[["iter", "STRING_DENS"]].dropna(how="any").values
    ax3.plot(vals[:, 0], vals[:, 1], **plot_kwargs)
    ax3.set_ylabel('Avg STRING Network Density')
    ax3.set_xlabel('Iteration')

    ax4 = fig.add_subplot(3,3,4, sharex=ax1)
    vals = stats[["iter", "r0"]].dropna(how="any").values
    plot_kwargs = dict(marker=',')
    ax4.plot(vals[:, 0], vals[:, 1], color='k', **plot_kwargs)
    vals = stats[["iter", "m0"]].dropna(how="any").values
    ax4.plot(vals[:, 0], vals[:, 1], color='r', **plot_kwargs)
    vals = stats[["iter", "n0"]].dropna(how="any").values
    ax4.plot(vals[:, 0], vals[:, 1], color='g', **plot_kwargs)
    ax4.set_xlabel('Iteration')

    text_kwargs = dict(transform=plt.gca().transAxes,
                       ha='center', va='bottom')
    ax4.text(0.20, 0.95, 'Row scaling', color='k', **text_kwargs)
    ax4.text(0.50, 0.95, 'Mot scaling', color='r', **text_kwargs)
    ax4.text(0.80, 0.95, 'Net scaling', color='g', **text_kwargs)
    
    ax2a = fig.add_subplot(3,3,5, sharex=ax1)
    vals = stats[["iter", "ROWS"]].dropna(how="any").values
    ax2a.plot(vals[:, 0], vals[:, 1], **plot_kwargs)
    ax2a.set_ylabel('Avg Num Genes')
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
    ax4.text(0.20, 0.95, 'Rows', color='k', **text_kwargs)
    ax4.text(0.80, 0.95, 'Cols', color='r', **text_kwargs)

    ax5 = fig.add_subplot(3,3,7)
    vals = stats[["iter", "N_MOVES"]].dropna(how="any").values
    ax5.plot(vals[:, 0], vals[:, 1], color='k', **plot_kwargs)
    vals = stats[["iter", "N_IMPROVEMENTS"]].dropna(how="any").values
    ax5.plot(vals[:, 0], vals[:, 1], color='r', **plot_kwargs)
    text_kwargs = dict(transform=plt.gca().transAxes,
                       ha='center', va='bottom')
    ax5.text(0.20, 0.95, 'Num Moves', color='k', **text_kwargs)
    ax5.text(0.70, 0.95, 'Num Improvements', color='r', **text_kwargs)
    ax5.set_xlabel('Iteration')
    ax5.set_ylabel('Count')

    ax5a = fig.add_subplot(3,3,8)
    resids = np.array( [c.resid for c in globals.clusters.values()] )
    resids = resids[ resids > 0.001 ]
    plt.hist( resids, 20 )
    ax5a.set_xlabel('Cluster Residuals')
    ax5a.set_ylabel('Count')

    ax6 = fig.add_subplot(3,3,9)
    pclusts = np.array( [c.meanp_meme for c in globals.clusters.values()] )
    pclusts = pclusts[ ~np.isnan(pclusts) & ~np.isinf(pclusts) ]
    plt.hist( pclusts, 20 )
    ax6.set_xlabel('Cluster Mean log10(P-value)s')
    ax6.set_ylabel('Count')

    plt.show()


if __name__ == '__main__':
    import init
    init.init( 'output/Hpy.pkl' )
    
    plot_stats()
