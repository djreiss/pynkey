## DONE: handle fasta files w/ multiple entries. Right now, concatenates them all into one seq.
## TODO: Correctly handle IUPAC
## DONE: Need to add remove-ATG step, and also running dust on sequences.

import os
import tempfile
import warnings

import pandas as pd
import numpy as np

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio import SeqIO

import params

print 'sequence'

def get_genome_seq( genome_seqs, scaffoldId ):
    sid = str(scaffoldId)
    return genome_seqs[sid].seq ## it is a SeqRecord - a dict with a string as the key
    # genome_seq_names = np.array( [ seq.name for seq in genome_seqs.values() ] )
    # ind = np.where( genome_seq_names == sid )[ 0 ]
    # if len(ind) > 0:
    #     return genome_seqs.values()[ind[0]].seq
    # else:
    #     return ''

def get_sequences( genes, anno=None, ##=globals()['anno'], 
                   genome_seqs=None, ##=globals()['genome_seqs'], 
                   op_table=None, ##=globals()['op_table'], 
                   op_shift=True, distance=params.distance_search, debug=False ):
    # if anno is None:
    #     anno = globals()['anno']
    # if genome_seqs is None:
    #     genome_seqs = globals()['genome_seqs']
    # if op_table is None:
    #     op_table = globals()['op_table']

    seqs = {}
    ind = 0
    for gene in genes:
        upstream_gene = gene
        if gene not in anno.index: ##anno_ind == 0 ##! any(anno["sysName"].data .== gene) 
            seqs[ gene ] = pd.DataFrame( {'gene':[gene], 'upstream_gene':[upstream_gene], 'seq':[''], 'strand':['']}, 
                                         index=['gene'] )
            continue
        gene_anno = anno.ix[gene, :]
        if op_shift and op_table.shape[0] > 0:
            strand = gene_anno['strand']
            upLabel = 'SysName1' if strand == '+' else 'SysName2'
            currLabel = 'SysName2' if strand == '+' else 'SysName1'
            row = np.where( op_table[currLabel] == gene )[0]
            while len(row) > 0 and op_table.ix[ row[0], "pOp" ] >= 0.8:
                upstream_gene = op_table.ix[ row, upLabel ].values[0]
                row = np.where( op_table[currLabel] == upstream_gene )[0]

        tmp_anno = anno.ix[upstream_gene, :]
        strnd = tmp_anno['strand']
        strt = tmp_anno['start']
        rng = strt + distance if strnd == '+' else strt - np.array( [distance[1], distance[0]] ) ## distance.reverse() actually reverses the array inplace.
        genome_seq = get_genome_seq( genome_seqs, tmp_anno['scaffoldId'] )
        if rng[1] > len(genome_seq):
            rng[1] = len(genome_seq)
        seq = genome_seq[ rng[0]:rng[1] ]
        #if debug: print gene + ' ' + upstream_gene + ' ' + strnd + ' ' + strt + ' ' + rng + ' ' + seq
        if strnd == "-":
            seq = seq.reverse_complement()
        #print gene
        seqs[ gene ] = pd.DataFrame( {'gene':[gene], 'upstream_gene':[upstream_gene], 'seq':[seq.tostring()], 
                                      'strand':[strnd]}, index=['gene'] )

    #print seqs.keys()
    print len(seqs); print(seqs[genes[0]])
    out = pd.concat(seqs, ignore_index=True)
    out.index = out.gene
    return out

def filter_sequences( seqs, distance=params.distance_search, remove_repeats=True, remove_atgs=True ):
     seqs = seqs[ seqs.seq != '' ] ## remove all empty sequences; they cause heartburn.

     if remove_repeats: ##&& len( grep( "NNNNNN", seqs ) ) <= 1
         ##if verbose: println( "Removing low-complexity regions from sequences.\n" )
         fname = tempfile.mktemp() 
         writeFasta( seqs, fname )

         tmp = os.popen( './progs/dust %s' % fname ).read() ## note dust has an optional numeric parameter -- for what?
         os.unlink( fname )
         fname = tempfile.mktemp()
         handle = open(fname, 'w')
         handle.write(tmp)
         handle.close()

         handle = open(fname, 'r')
         seqs_new = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
         handle.close()
         os.unlink( fname )

         if seqs.shape[0] != len(seqs_new):
             warnings.warn("Remove low complexity failed - skipping!" )
         else:
             out = seqs.copy()
             for gene in seqs.gene:
                 newseq = seqs_new[gene].seq.tostring()
                 out.seq[gene] = newseq
             seqs = out

     if remove_atgs and np.any( distance < 0 ): 
         nnnn = 'NNNN'
         for i in seqs.index:
             ss = seqs.seq[i]
             ##ss[ (distance[1]+0):(distance[1]+4) ] = nnnn  ## Mask out ATGs (why 4 instead of 3?)
             if seqs.strand[i] == '-':
                 ss = ss[:(distance[1]-1)] + nnnn + ss[(distance[1]+3):]
             elif seqs.strand[i] == '+':
                 ss = ss[:(len(ss)-distance[1]-1)] + nnnn + ss[(len(ss)-distance[1]+3):]
             ##print i; print seqs.strand.values[i]; print seqs.seq.values[i]; print ss
             seqs.seq[i] = ss 
     return seqs

## Create DataFrame for multiple sequences; same format as get_sequences(bicluster)
## TODO: possibly store strand and upstream_gene info in fasta file for fetching later
def readFastaDNA( fname ):
    handle = open(fname, 'r')
    seqs = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
    handle.close()

    out = [ pd.DataFrame( {'gene':[s.id], 'upstream_gene':[''], 'seq':[s.seq.tostring()], 'strand':['']}, 
                           index=['gene'] ) for s in seqs.values() ]
    out = pd.concat(out, ignore_index=True)
    out.index = out.gene
    return out

def writeFasta( seqs, fname ): ## Assumes seqs in DataFrame format of get_sequences()
    tmp = [ SeqRecord(Seq(seqs.seq.values[i], IUPAC.ExtendedIUPACDNA()), id=seqs.gene.values[i], 
                      name=seqs.gene.values[i], description=seqs.gene.values[i]) for i in range(len(seqs)) ]
    handle = open(fname, 'w')
    SeqIO.write(tmp, handle, "fasta")
    handle.close()

def revComp( seq ):
    Seq.reverse_complement(seq)

DNA_letters = np.array( [ 'G', 'A', 'T', 'C' ] )

def generateAllKmers( length ):
    from sklearn.utils.extmath import cartesian
    all_combos = pd.DataFrame( cartesian( [DNA_letters] * length ) )
    all_combos = all_combos.apply( lambda x: ''.join( x.tolist() ), axis=1 ) ## sweet!
    return all_combos.values

# This is sweet!!!
def getBgCounts( seqs, order=3, include_revComp=True ):
    from sklearn.utils.extmath import cartesian
    d = {}
    for ord in range(order+1): ## get counts for 1,2,3,...order
        all_combos = generateAllKmers( ord+1 )
        for ss in all_combos:
            if ss not in d.keys():
                d[ss] = 0
            for i in range( len(seqs) ):
                if type(seqs) == np.ndarray:
                    seq = Seq(seqs[i], IUPAC.ExtendedIUPACDNA())  ## seqs is a column.values from a DataFrame 
                elif type(seqs) == dict: 
                    seq = seqs.values()[i].seq ## otherwise assumed to be a dict of Bio.Seqs (e.g. genome_seqs)
                d[ss] += seq.count(ss)
                d[ss] += seq.reverse_complement().count(ss)
    return d

## Convert counts dictionary into frequencies; divide counts for each k-mer by the total for that given k.
def getBgFreqs( bgCounts ):
    k = np.array(bgCounts.keys())
    nc = np.array([ len(i) for i in k ])
    d = {}
    for i in range( nc.max() ):
        ks = k[ nc == i+1 ]
        tot = sum( np.array([bgCounts[j] for j in ks]) )
        for j in ks:
            d[j] = float(bgCounts[j]) / float(tot)
    return d

# function read_iupac()
#     fname = "./IUPAC-dna.txt"
#     str = open( fname )
#     lines = split( readall( str ), '\n' )
#     close( str )
#     lines = [ replace(l," ","") for l=lines]
#     lines = [ replace(l,r"\(.*\)","") for l=lines ]
#     lines = lines[ lines .!= "" ]
#     d = Dict()
#     for l in lines
#         if any(DNA_letters .== l[1]) d[l[1]] = l[1]; 
#         elseif l[1] == 'U' d[l[1]] = l[1]; 
#         else d[l[1]] = l[3:end]; end
#     end
#     d
# end

# const iupac_dict = read_iupac();
