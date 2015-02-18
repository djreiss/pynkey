## TODO: Try using biopython for sequences
## TODO: Load and use other networks (with their respective scores)
import os
import gzip

import numpy as np
from numpy import nan as NA
from numpy import random as rand
import pandas as pd
import networkx as nx

import sequence as seq
import params
from Bicluster import bicluster
import globals as glb

print 'importing init'

IS_INITED = False

def init( from_pickle_file=None ):
    global IS_INITED
    #global ratios, genome_seqs, anno, op_table, string_net, allSeqs_fname, all_bgFreqs, all_genes
    #global pynkey_code, clusters, endTime

    if from_pickle_file is not None:
        from funcs import load_checkpoint
        import logging
        logging.info( 'Loading checkpoint file: %s' % from_pickle_file )
        load_checkpoint( from_pickle_file )
        print str(glb.startTime)

    else:

        glb.ratios, glb.genome_seqs, glb.anno, glb.op_table, glb.string_net, glb.allSeqs_fname, \
            glb.all_bgFreqs, glb.all_genes = \
            pynkey_init(params.organism, params.k_clust, params.ratios_file)

        if params.k_clust <= 0: ## default for organisms not Eco/Hpy/Mpn/Hal/Sce based on total # of genes
            params.k_clust = int( len(glb.all_genes) / params.avg_genes_per_cluster * params.avg_clusters_per_gene )
            print params.k_clust

        # Save all pynkey code for safe keeping
        glb.pynkey_code = {}
        try:
            glb.pynkey_code = load_pynkey_code()
        except:
            print 'Must have "import"ed this module rather than run from main.py.'
        
        glb.clusters = init_biclusters( glb.ratios, params.k_clust, 'kmeans+random' );

        ##counts_g = bicluster.get_all_cluster_row_counts( clusters, all_genes )

    import datetime
    print 'DONE WITH INITIALIZATION!'
    glb.endTime = datetime.datetime.now()
    print str(glb.endTime)
    print str(glb.endTime - glb.startTime) + ' seconds since initialization'
    IS_INITED = True

def pynkey_init(organism, k_clust, ratios_file):
    print organism

    ratios = load_ratios(ratios_file)

    genome_seqs = load_genome(organism)

    anno = load_annos(organism)

    string_net = load_string_net(organism)

    op_table = load_op_table(organism)

    ## build up a list of all genes in the expr. data + in the annotations + in the string network
    all_genes = ratios.index.values ## same as rownames
    if params.all_genes_option == 'all': ## append all gene names found in other data types
        all_genes = np.unique( np.append( all_genes, anno.index.values ) )
        all_genes = np.unique( np.append( all_genes, [string_net.protein1.values, string_net.protein2.values] ) )
        all_genes = np.unique( np.append( all_genes, [op_table.SysName1.values, op_table.SysName2.values] ) )
    else: ## only use gene names found in expression data; remove it from the networks at least
        ##string_net = string_net[ np.in1d( string_net.protein1, all_genes ) & \
        ##                         np.in1d( string_net.protein2, all_genes ) ]
        string_net = string_net.subgraph(all_genes)

    all_genes = np.sort(all_genes)
    all_genes = all_genes[ all_genes != NA ] ## for some reason some are floats and they are printing as nan's but this line
    if type(all_genes[0]) == float:          ## does not remove them so we do this line instead
        all_genes = all_genes[1:]
    print all_genes.shape

    gene_regex = get_regex( all_genes )
    print gene_regex
    
    ## Get all upstream seqs and their bgFreqs
    ## DONE: Need to add 0's for k-mers that are NOT in the sequences.
    ## TODO: Need to include vague IUPAC symbols better
    print 'Getting upstream sequences (search)'
    all_seqs = seq.get_sequences(all_genes, anno, genome_seqs, op_table, True, params.distance_search, False) 
    print 'Filtering upstream sequences'
    all_seqs = seq.filter_sequences(all_seqs, params.distance_search)

    print 'Getting upstream sequences (scan)'
    all_seqs_scan = seq.get_sequences(all_genes, anno, genome_seqs, op_table, True, params.distance_scan, False)
    print 'Filtering upstream sequences'
    all_seqs_scan = seq.filter_sequences(all_seqs_scan, params.distance_scan)

    allSeqs_fname = './%s/allSeqs.fst' % params.organism
    seq.writeFasta( all_seqs_scan, allSeqs_fname ) ## NOTE all_seqs_scan are not used from here on, just the fasta file
    ## NOTE: can use the fasta file as an offline database:
    ## test=SeqIO.index("Hpy/allSeqs.fst", "fasta") ## stores a dict with links to the locations in the file
    ## test['HP0001'] ## works!!
    
    # ##bgCounts = getBgCounts( [genome_seqs,revComp(genome_seqs)], [0:5], true );
    # ##bgFreqs = getBgFreqs( bgCounts ); ## TODO: bgFreqs are currently not used in MEME-ing OR MAST-ing.

    print 'Getting trinucleotide freqs from upstream sequences'
    all_bgCounts = seq.getBgCounts( all_seqs.seq.values )
    all_bgFreqs = seq.getBgFreqs( all_bgCounts );  ## TODO: bgFreqs are currently not used in MEME-ing OR MAST-ing.

    # save_jld( "./$organism/data.jldz", (organism, k_clust, ratios, genome_seqs, anno, op_table, string_net, 
    #                                       allSeqs_fname, all_bgFreqs, all_genes) ) ##, all_rows
    
    # (ratios, genome_seqs, anno, op_table, string_net, ##all_seqs, all_seqs_scan, 
    #  allSeqs_fname, all_bgFreqs, all_genes) ##, all_rows) ##all_bgCounts, 
    return (ratios, genome_seqs, anno, op_table, string_net, allSeqs_fname, all_bgFreqs, all_genes)

def load_ratios(rats_file):
    ## TODO: make it work for any organism; download the sequence/anno data using HTTP package; 
    ## DONE: try on bigger data; multiple chrome's
    print rats_file
    x = pd.read_table(rats_file, compression='gzip', index_col=0) ## note x.ix[:,:5].head() prints first 5 cols
    print x.shape 

## Examples for Eco:
## First get the taxinomy codes:
## wget -O taxdump.tar.gz 'ftp://ftp.ncbi.nih.gov/pub/taxonomy/taxdump.tar.gz'
## tar -xvzf taxdump.tar.gz names.dmp
## wget -O genomeInfo.511145.tsv 'http://www.microbesonline.org/cgi-bin/genomeInfo.cgi?tId=511145;export=tab'
## wget -O genome.511145.txt 'http://microbesonline.org/cgi-bin/genomeInfo.cgi?tId=511145;export=genome'
## wget -O STRING_511145.tsv.gz 'http://baliga.systemsbiology.net/data/STRING/511145_STRING.tsv.gz'; gunzip STRING_511145.tsv.gz
## NOTE: extract first 200 columns from a file using cut:
##    gunzip -c TanayMSB_GEData.txt.gz | cut -f1-200 | more

    ##means=x.apply(mean,1)
    ##sds=x.apply(std,1)

    ## Note: ratios.as_matrix() converts to a numpy array
    ## Although ratios.values does the same thing ;)
    ## Note e.g. ratios.mean(1) is about 5.5x slower than ratios.values.mean(1)

    good_rows = x.isnull().apply(np.sum, axis=1) < x.shape[1]/3
    x = x[good_rows]
    print x.shape

    ##for i=1:size(x,1) x[i,:] = stdize_vector(x[i,:]); end ##(x[i,:]-nanmean(x[i,:]))/nansd(x[i,:]); end
    f = lambda x: (x-x.dropna().mean()) / x.dropna().std()
    x = x.apply(f, axis=1)
    print x.ix[:, :6].head()
    return x

def load_genome(organism):
    from Bio import SeqIO
    ## see http://biopython.org/wiki/SeqIO#Sequence_Input
    ## DONE: allow for multiple genome seqs (e.g. Halo, Sce)
    org_files = np.array( os.listdir('./' + organism + '/') )
    genome_file = './' + organism + '/' + org_files[ np.array( [f.startswith('genome.') for f in org_files] ) ][ 0 ]
    print genome_file

    handle = gzip.open(genome_file)
    genome_seqs = SeqIO.to_dict(SeqIO.parse(handle, 'fasta')) ## allows >1 sequence
    handle.close()
    ##genome_seqs = SeqIO.read(gzip.open(genome_file), "fasta") ## this is a SeqRecord (only one sequence)
    ## get the biosequence via: genome_seqs.values()[0].seq

    ## NOTE can avoid reading in sequence into memory! 
    ## See http://biopython.org/DIST/docs/tutorial/Tutorial.html#sec%3ASeqIO-index
    ## genome_seqs = SeqIO.index(genome_file, 'fasta') ## probably need to ungz the file first.

    print genome_seqs
    print len(genome_seqs)

    ## TODO: can use carray to compress the genome seq in memory:
    ##b = carray.carray(np.array(list(str(genome_seqs.values()[0]))))
    return genome_seqs

def load_annos(organism):
    import string ## for replace

    ## Load the gene annotations
    org_files = np.array( os.listdir('./' + organism + '/') )
    genomeInfo_file = './' + organism + '/' + \
        org_files[ np.array( [f.startswith('genomeInfo.') for f in org_files] ) ][ 0 ]

    print genomeInfo_file
    ## note x.ix[:,:5].head() prints first 5 cols
    x = pd.read_table(genomeInfo_file, compression='gzip', index_col='sysName') 
    x['desc'].values[np.where(np.array([type(s)==str for s in x['desc'].values])==False)] = '' ## replace nan's w ''
    x['desc']=np.array([string.replace(s,' (NCBI ptt file)','') for s in x['desc'].values]) ## makes it prettier
    return x

## TBD: save string_net as a sparse data frame (i.e., matrix) for speed?
## see http://pandas.pydata.org/pandas-docs/dev/sparse.html
def load_string_net(organism):
## Load the string network
    org_files = np.array( os.listdir('./' + organism + '/') )
    string_file = './' + organism + '/' + org_files[ np.array( [f.startswith('STRING') for f in org_files] ) ][ 0 ]

    print string_file
    ##string_net = pd.read_table(string_file, compression='gzip', names=['protein1','protein2','weight'], 
    ##                           index_col=False) 
    ##string_net.set_index( 'protein1', drop=False, inplace=True )

    ## Symmetrize it? -- no, seems to already be done.
#   tmp = pd.concat([string_net.ix[:,1], string_net.ix[:,0], string_net.ix[:,2]], axis=1, ignore_index=True)
#   tmp.columns = tmp.columns[[1,0,2]] ## reorder the columns to the same as string_net
#   string_net = pd.concat( [string_net, tmp], ignore_index=True )

    ##print string_net.shape

    string_net = nx.read_weighted_edgelist( string_file )

    return string_net

def load_op_table(organism):
## Now, require the operons table!
## DONE? (need to verify): Allow for no operons table (e.g. yeast)
## TODO: reindex op_table by SysName1 -- choosing via index is about 4x faster than [np.where(...)]
    org_files = np.array( os.listdir('./' + organism + '/') )
    op_table = pd.DataFrame()
    try:
        op_file = './' + organism + '/' + org_files[ np.array( [f.startswith('microbesonline_operons_') 
                                                                for f in org_files] ) ][ 0 ]
        print op_file
        op_table = pd.read_table(op_file, compression='gzip') 
        op_table = op_table[ [ 'SysName1', 'SysName2', 'bOp', 'pOp' ] ]
        print op_table.shape
    except:
        print 'No operons file'
        
    return op_table

## Load the pynkey code as text so it can be stored in the run results for rerunning
## from here: http://code.activestate.com/recipes/496889-dynamically-determine-execution-path-of-a-file/
import os, sys, inspect
def load_pynkey_code():
    def execution_path(filename):
        return os.path.join(os.path.dirname(inspect.getfile(sys._getframe(1))), filename)
    exec_path = execution_path('') 
    print exec_path
    files = np.array( os.listdir(exec_path) )
    files = files[ np.array( [f.endswith('.py') for f in files] ) ]
    files = files[ np.array( [not f.startswith('.') for f in files] ) ]
    print files
    code = {}
    for f in files:
        fo = open( exec_path + '/' + f, 'r' )
        lines = fo.readlines()
        fo.close()
        code[f] = lines
    return code

def get_regex( strings, min_ignore=2 ):
    nchar = np.array( [ len(i) for i in strings ] )
    out = ''
    for i in xrange( nchar.max() ):
        d = {}
        for str in strings:
            if i >= len(str):
                continue
            c = str[i]
 	    if c in d:
                d[c] = d[c] + 1
            else:
                d[c] = 1
        if len(d) == 0:
            break
        tmp = ''
        for c in d:
            if d[c] > min_ignore: ## ignore values with <=2 occurences
                tmp = tmp + c ##''.join([tmp, c]) ##'{0}{1}'.format(tmp, c)
        if len(tmp) > 1:          ## add brackets around it if >1 different characters
            tmp = '[' + tmp + ']' ##''.join(['[', tmp, ']']) ##'[{0}]'.format(tmp)
        if sum(d.values()) < len(strings) * 0.95:  ## add '?' after it if infrequent
            tmp = tmp + '?' ##''.join([tmp, '?']) ##'{0}?'.format(tmp)
        if tmp != '':             ## append it
            out = out + tmp ##''.join([out, tmp]) ##'{0}{1}'.format(out, tmp)
    return out

def init_biclusters( ratios, k_clust, method='kmeans' ):
    import scipy.cluster.vq as clust

## Init via kmeans -- TODO: other methods (including random)
    clusters = {}

    print method
    if method == 'kmeans' or method == 'kmeans+random':
        nans = np.sum(np.isnan(ratios.values))
        x = ratios
        if nans > 0:
            x = ratios.copy()
            ## can kemans2 handle na's?
            ##x.fillna(method=None, value=rand.normal(scale=0.1, size=nans), axis=1, inplace=True)
            x.values[ x.isnull().values ] = rand.normal(scale=0.1, size=nans)

        _, km1 = clust.kmeans2( x.values, k_clust, iter=20, minit='random' )

         ## seed each bicluster with rows=output from kmeans and cols=random (1/2 of all cols)
        for k in xrange(k_clust):
            rows = ratios.index.values[ np.where( km1 == k ) ]
            if method == 'kmeans': 
                if len(rows) == 0:
                    rows = ratios.index.values[ rand.choice(ratios.shape[0], 10, replace=False) ]
                clusters[k] = bicluster( k, rows, x )
            elif method == 'kmeans+random':          ## add some random rows to each bicluster
                new_rows = ratios.index.values[ rand.choice(ratios.shape[0], max(len(rows),10), replace=False) ]
                clusters[k] = bicluster( k, np.unique( np.concatenate( (rows, new_rows), 0 ) ), x )
    elif method == 'random':
        for k in xrange(k_clust):
            rows = ratios.index.values[ rand.choice(ratios.shape[0], 20, replace=False) ]
            clusters[k] = bicluster( k, np.unique( rows ), x )
    return clusters

