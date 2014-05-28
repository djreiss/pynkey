## TODO: handle fasta files w/ multiple entries. Right now, concatenates them all into one seq.
## TODO: Correctly handle IUPAC
## TODO: Need to add remove-ATG step, and also running dust on sequences.

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

def get_genome_seq( genome_seqs, scaffoldId ):
    sid = str(scaffoldId)
    genome_seq_names = np.array( [ seq.name for seq in genome_seqs.values() ] )
    ind = np.where( genome_seq_names == sid )[ 0 ]
    if len(ind) > 0:
        return genome_seqs.values()[ind[0]].seq
    else:
        return ''

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

def seqs_to_seqRecord( seqs, names ):
    return [ SeqRecord(Seq(seqs[i], IUPAC.ExtendedIUPACDNA), id=names[i]) for i in range(len(seqs)) ]

def filter_sequences( seqs, distance=params.distance_search, remove_repeats=True, remove_atgs=True ):
     seqs = seqs[ seqs.seq != '' ] ## remove all empty sequences; they cause heartburn.

     if remove_repeats: ##&& len( grep( "NNNNNN", seqs ) ) <= 1
         ##if verbose: println( "Removing low-complexity regions from sequences.\n" )
         ## Need to convert array of strings back to sequence record for output via biopython:
         tmp = [ SeqRecord(Seq(seqs.seq.values[i], IUPAC.ExtendedIUPACDNA), id=seqs.gene.values[i], 
                           name=seqs.gene.values[i], description=seqs.gene.values[i]) for i in range(len(seqs)) ]
         fname = tempfile.mktemp() 
         handle = open(fname, 'w')
         SeqIO.write(tmp, handle, "fasta")
         handle.close()

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

# ## Create DataFrame for multiple sequences; same format as get_sequences(bicluster)
# function readFastaDNA( fname ) 
#     str = open( fname )
#     seq::ASCIIString = readall( str );
#     close( str );
#     seqs = split( seq, '>' );
#     first_return = [search(seqs[i], '\n') for i=2:length(seqs)] ## first seq is ""
#     seq_names = convert( Vector{ASCIIString}, [seqs[i][1:(first_return[i-1]-1)] for i=2:length(seqs)] ); 
#     seqs = [seqs[i][(first_return[i-1]+1):end] for i=2:length(seqs)];
#     seqs = convert( Vector{ASCIIString}, [uppercase( replace( replace( seqs[i], r">.*\n", "" ), '\n', "" ) ) for i=1:length(seqs)] );
#     out = DataFrame( {"gene"=>seq_names, "upstream_gene"=>fill("", length(seqs)), "seq"=>seqs} );
#     out
# end

# function writeFasta( seqs, fname ) ## Assumes seqs in DataFrame format of get_sequences()
#     str = open( fname, "w" )
#     for i=1:size(seqs,1)
#         gene = seqs["gene"].data[i]
#         seq = seqs["seq"].data[i]
#         write( str, ">$gene\n" )
#         write( str, "$seq\n" )
#     end
#     close( str )
# end

# function revComp( seq::ASCIIString )
#     out = Array( Uint8, length( seq ) ) ## uninitialized is faster?
#     j = length( seq ) + 1
#     seq = uppercase( seq )
#     for i = 1:length( seq )
#         c = seq[ i ]
#         out[ j -= 1 ] = c == 'G' ? 'C' : ( c == 'C' ? 'G' : ( c == 'T' ? 'A' : ( c == 'A' ? 'T' : c ) ) )
#     end
#     ASCIIString( out );
# end

# const DNA_letters = [ 'G', 'A', 'T', 'C' ];
# const DNA_letter_lookup = {'G'=>1, 'A'=>2, 'T'=>3, 'C'=>4 };

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

# # function collect_add_dicts( a, b ) ## Add up all elements into output dict
# #     for (k,v) in b a[k] = v + b[k]; end
# #     a
# # end

# #@profile begin
# ## Note this will not have an entry for subseqs that don't exist in the sequence
# ## Only does it for a SINGLE "order" value and single strand - see the default getBgCounts() to see standard usage.
# ## Use in conjunction with generate_all_kmers in order to get seqs w/ count=0
# ## Taken from my seqTest.jl
# function getBgCounts( seqs::Array{ASCIIString,1}, order::Array{Int64,1}=[0:5], verbose::Bool=false ) 
#     ss::ASCIIString = "";
#     d = Dict{ASCIIString,Int64}() ## Specifying the types speeds it up about 35%
#     for seq=[seqs,[revComp(s) for s=seqs]] ##seqs
#         seq = uppercase( seq )
#         for j=order
#             if verbose println(j); end
#             for i=1:length(seq)-j ## for loop is about 5% faster than while loop
# 	        ss = seq[i:(i+j)]
# 	        d[ss] = get(d,ss,0) + 1
#             end
#         end
#     end
#     d
# end

# ##getBgCounts( seqs::Array{ASCIIString,1}, order::Int64 ) = 
# ##             getBgCounts( [seqs,[revComp(s) for s=seqs]], [order] );
# # getBgCounts( seqs::Array{ASCIIString,1}, order::Array{Int64,1} ) = 
# #                getBgCounts( [seqs,[revComp(s) for s=seqs]], order ); ##, false );
# # getBgCounts( seqs::Array{ASCIIString,1} ) = 
# #                getBgCounts( [seqs,[revComp(s) for s=seqs]] );
# # getBgCounts( seq::ASCIIString ) = getBgCounts( [seq] );

# #bgCounts = getBgCounts( genome_seqs );
# #end # profile

# ## Convert counts dictionary into frequencies; divide counts for each k-mer by the total for that given k.
# function getBgFreqs( bgCounts::Dict{ASCIIString,Int64} )
#     k = collect(keys( bgCounts ))
#     nc = [ length(k[i]) for i=1:length(k) ]
#     d = Dict{ASCIIString,Float64}()
#     for i=1:maximum(nc)
#         ks = k[find(nc.==i)]
#         tot = sum( [bgCounts[j] for j=ks] )
#         for j=ks d[j] = bgCounts[j] / tot; end
#     end
#     d
# end

# ## Couldn't do this in R! Generate all sequences with given length.
# ## Allow for up to n_ns N's as well.
# function generate_all_kmers( len::Int64, n_ns::Int64=0, letters::Array{Char,1}=DNA_letters )
#     function append_nucs( d::Array{ASCIIString}, letters::Array{Char,1} )
#         out = Array(ASCIIString, length(d)*length(letters))
#         ind::Int64 = 0
#         for i=d for j=letters out[ind+=1] = strcat(i,j); end; end ## note can also use i*j
#         out
#     end

#     if n_ns > 0 letters = [ letters, 'N' ]; end ## Add N to letters; just need 1 copy.
#     d = Array(ASCIIString, length(letters))
#     i=0; for j=letters d[i+=1] = "$j"; end
#     for i = 1:len-1 d = append_nucs( d, letters ); end
#     if n_ns > 0   ## Get rid of k-mers with > n_ns Ns.
#         n_count = [ sum(chars(d[i]) .== 'N') for i=1:length(d) ]
#         d = d[ find( n_count .<= n_ns ) ]
#     end
        
#     d
# end

# #generate_all_kmers( len::Int64 ) = generate_all_kmers( len, 0, DNA_letters )
# #generate_all_kmers( len::Int64, n_ns::Int64 ) = generate_all_kmers( len, n_ns, DNA_letters )

# if false
#     seq = "GATCATGCATGTATGCTACGTGCGCGGGTACGTATATGATGCTATTATCGTAGCTACGTAGCTAGCTAGCTACAGTCGATCGATTGAC"
#     seq2bit = dna2seq(seq)
#     d=Dict{BitArray{1},Int64}()
#     d[[seq2bit[1:5].b1,seq2bit[1:5].b2]] = 1 ## can count up k-mers same as with ASCIIString

#     ## file-based - can handle IUPAC codes
#     seqnt = nt(seq)
#     Astream=open("test.mmap", "a+")
#     A=mmap_array(Nucleotide, size(seqnt), Atream)
#     A[1:length(seqnt)] = seqnt
#     close(Astream)

#     ## try this! file-based, don't even have to read it in!
#     len = filesize("Hpy/genome.85962.txt")
#     seqStream = open("Hpy/genome.85962.txt", "r")
#     seq = mmap_array(Nucleotide, (len,), seqStream)
#     println(convert(ASCIIString,seq[500:505]))
# end
