## DONE: allow parallelization by adding seqs to the bicluster struct, sending biclusters (and allSeqs_fname) 
##    through to all child nodes. Need to make a function that takes a single argument (Dict?) and use pmap().
## DONE: Add n_motifs to the schedule so it ups to 2 after n_iters*3/4 iterations -- NOTE needs to get passed to
##    child processes
## DONE: add min/max_width and search/scan_distance paramters as globals (e.g. to change for yeast)
## DONE: add filter_sequences function to run dust and remove ATGs

import tempfile
import os
import warnings
import sys

import numpy as np
import pandas as pd

import sequence as seq

print 'importing meme'

def get_n_motifs( iter, n_iters ):
    n_motifs = 1 if ( iter < n_iters * 3 / 4 ) else 2 ## n_motifs is 1 in early iterations; 2 in later.
    return n_motifs

## DONE: only send (and receive) the meme_out and mast_out data (not entire bicluster!) -- a bit faster..?
# function re_meme_all_biclusters_parallel( clusters::Dict{Int64,bicluster}, force::Bool=false, verbose::Bool=false )
#     global k_clust, allSeqs_fname, iter, n_iters, distance_search
#     data::Array{Any,1} = []
#     n_motifs = get_n_motifs( iter, n_iters )
#     r_rownames = rownames(ratios)
#     for k in 1:k_clust
#         b = clusters[k]
#         if ! force && ! b.changed[1] continue; end ## If rows not changed, skip it
#         seqs = get_sequences( r_rownames[b.rows] ) 
#         seqs = filter_sequences( seqs, distance_search )
#         dat = { "k" => b.k, "seqs" => seqs, "allSeqs_fname" => allSeqs_fname, "n_motifs" => n_motifs, "verbose" => verbose };
#         push!( data, dat )
#     end
#     memeOut = pmap( re_meme_bicluster, data ) ## returns an Array{Any,1}
#     for i in 1:length(memeOut)
#         if memeOut[i][2] == nothing continue; end
#         k = memeOut[i][1]
#         clusters[k].meme_out = memeOut[i][2]
#         clusters[k].mast_out = memeOut[i][3]
#     end
#     clusters
# end

# ## for parallelization from re_meme_all_clusters_parallel() :
# function re_meme_bicluster( x::Dict{Any,Any} ) 
#     if x["seqs"] == nothing 
#         k = x["k"] ## shouldn't happen!
#         warn( "NO GENES! $k" )        
#         return( (k, nothing, nothing) )
#     end 
#     re_meme_bicluster(x["k"], x["seqs"], x["n_motifs"], x["allSeqs_fname"], x["verbose"])
# end

# function re_meme_all_biclusters( clusters::Dict{Int64,bicluster}, force::Bool=false )
#     global k_clust, iter, n_iters
#     n_motifs = get_n_motifs( iter, n_iters )
#     for k in 1:k_clust
#         b = clusters[k]
#         if ! force && ! b.changed[1] continue; end ## If rows not changed, skip it
#         b = re_meme_bicluster( b, n_motifs )
#         clusters[k] = b
#     end
#     clusters
# end

## Called from bicluster.re_meme()
def re_meme_bicluster( k, seqs, n_motifs, allSeqs_fname, motif_width_range, verbose=False ):
    #gibbs_out = gibbs_site_sampler(seqs["seq"].data)     ## run gibbs sampler on most "flagellar-enriched" cluster
    #gibbs_out2 = gibbs_site_sampler(seqs, gibbs_out["pssm"])

    meme_out = []
    mast_out = pd.DataFrame()
    try:
        ns = seqs.shape[0] 
        if verbose:
            print 'IN HERE: MOT', n_motifs, k, ns
        meme_out = do_meme( seqs, motif_width_range, n_motifs, verbose )
        if meme_out != "": # NOTE for now if meme fails (most often b/c too few seqs) just keep the previous run
            try:
                mast_out = do_mast(meme_out, allSeqs_fname, False, verbose)
            except:
                warnings.warn( 'ERROR RUNNING MAST FOR BICLUSTER %d' % k )
                print sys.exc_info()[0]
        else:
            warnings.warn( 'TOO FEW SEQUENCES TO MEME FOR BICLUSTER %d' %k )
    except:
        warnings.warn( "TOO FEW SEQUENCES TO MEME FOR BICLUSTER %d" % k )
        print sys.exc_info()[0]
    return (k, meme_out, mast_out)

## This is parallelizable because it sends the seqs to be searched to each node
def do_meme(seqs, motif_width_range, n_motifs=2, verbose=False):
    seqs = seqs[ seqs.seq != '' ]
    seqs = seqs[ ~seqs.seq.duplicated().values ] ## note ~ is boolean not
    if seqs.shape[0]  <= 2:
        return ''

    seqs_fname = tempfile.mktemp() 
    seq.writeFasta( seqs, seqs_fname )
    
    cmd = './progs/meme %s -time 600 -dna -revcomp -maxsize 9999999 -nmotifs %d -evt 999999 -minw %d -maxw %d -mod zoops -nostatus -text' % (seqs_fname, n_motifs, motif_width_range[0], motif_width_range[1])
    if verbose:
        print cmd

    memeOut = os.popen( cmd ).read()
    os.unlink( seqs_fname )
    #if verbose: 
    #    print memeOut.split( '\n')
    return memeOut
                        
import StringIO
import re

## this IS parallelizable b/c it just uses the allSeqs_fname file for reading, which was already created
def do_mast(memeOut, allSeqs_fname, get_allHitsTab=False, verbose=True):
    memeOutFname = tempfile.mktemp()
    handle = open(memeOutFname, 'w')
    handle.write(memeOut)
    handle.close()

    ## Not used right now but possibly useful later for plotting motif positions
    if get_allHitsTab: ## First get the table of sites (multiple hits per sequence) -- using "-hits_list"
        cmd = './progs/mast %s -d %s -nostatus -stdout -text -brief -ev 999999 -mev 999999 -mt 1.0 -seqp -remcorr -hit_list' % (memeOutFname, allSeqs_fname)
        if verbose:
            print cmd
        mastOut = os.popen( cmd ).read()
        header = mastOut.split('\n')[1][2:].split(' ')
        tmpHitsTab = pd.read_csv( StringIO.StringIO(mastOut), sep=' ', names=header, skipinitialspace=True, 
                                  skiprows=2, skipfooter=1 )

    ## This gets the table of sequence e-values (one value per sequence) -- need to parse it out
    ## Table starts 2 lines below 
    ## SEQUENCE NAME                      DESCRIPTION                   E-VALUE  LENGTH
    ## and last line is 2 lines above 
    ## ********************************************************************************
    cmd = './progs/mast %s -d %s -nostatus -stdout -text -brief -ev 999999 -mev 999999 -mt 0.99 -seqp -remcorr' % (memeOutFname, allSeqs_fname) #-b
    if verbose: 
        print cmd

    mastOut = os.popen( cmd ).read()
    os.unlink( memeOutFname )

    ## We need P-values, not E-values, so we turned off the '-b' flag and need to parse the 3rd section (here is the 2nd section)
    ## We can get the e-vals from the 3rd section too, so let's use that
    tmp = np.array(mastOut.split('\n'))
    startInd = np.where( tmp == "SEQUENCE NAME                      DESCRIPTION                   E-VALUE  LENGTH" )[ 0 ][ 0 ] + 2
    endInd = np.where( tmp == "********************************************************************************" )[ 0]
    endInd = endInd[ endInd > startInd ][ 0 ] - 2
    mo = '\n'.join( tmp[ startInd:endInd ] )
    mastOutTab = pd.read_csv( StringIO.StringIO(mo), sep=' ', skipinitialspace=True, 
                              names=['Gene', 'E-value', 'Length'] )

    ## OK, now parse the 3rd section
    lines = np.where( np.array( [ i.find('COMBINED P-VALUE') for i in tmp ] ) >= 0 )[ 0 ]
    genes = tmp[ lines-2 ]
    reg = re.compile('COMBINED P-VALUE = (.+)(\s+E-VALUE = (.*))')
    ##p_values = np.array( [ reg.search( str(i) ).group(1,3) for i in tmp[lines] ] )
    p_values = np.array( [ float(reg.search( str(i) ).group(1).strip()) for i in tmp[lines] ] )
    e_values = np.array( [ float(reg.search( str(i) ).group(3).strip()) for i in tmp[lines] ] )
    mastOutTab = pd.DataFrame( {'Gene':genes, 'P-value':p_values, 'E-value':e_values} )
    
    if get_allHitsTab:
        return (mastOutTab, mastAllHitsTab)

    return mastOutTab 

