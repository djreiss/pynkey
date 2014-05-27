## TODO: Try using biopython for sequences
## TODO: Load and use other networks (with their respective scores)
from numpy import nan as NA
import numpy as np
import pandas as pd
import os
from Bio import SeqIO
import gzip

print 'init'

def pynkey_init(organism, k_clust, ratios_file):
    #global ratios_file, ratios, distance_scan, distance_search
    print organism

    ## TODO: shouldn't store ratios, k_clust in this file; this is just organism/genome data
    # if filesize("./$organism/data.jldz") > 0 ## 0 if not exist
    #     warn( "Loading organism data from ./$organism/data.jldz" )
    #     (organism, k_clust, ratios, genome_seqs, anno, op_table, string_net, 
    #      allSeqs_fname, all_bgFreqs, all_genes) = load_jld("./$organism/data.jldz"); ##, all_rows

    #     ## Make sure we have written out the sequences file for MAST-ing
    #     if filesize(allSeqs_fname) <= 0 
    #         all_seqs_scan = get_sequences(anno["sysName"].data,anno,genome_seqs,true,op_table,distance_scan,false); 
    #         all_seqs_scan = all_seqs_scan[ find(all_seqs_scan[:,1].!=""), : ]
    #         all_seqs_scan = filter_sequences( all_seqs_scan, distance_scan )

    #         writeFasta( all_seqs_scan, allSeqs_fname )
    #     end
    #     return( (ratios, genome_seqs, anno, op_table, string_net, ##all_seqs, all_seqs_scan, 
    #              allSeqs_fname, all_bgFreqs, all_genes) ) ##, all_rows) ##all_bgCounts, 
    # end 

    ratios = load_ratios(ratios_file)

    genome_seqs = load_genome(organism)

    anno = load_annos(organism)

    string_net = load_string_net(organism)

    op_table = load_op_table(organism)

    ## build up a list of all genes in the expr. data + in the annotations + in the string network
    all_genes = ratios.index.values ## same as rownames
    all_genes = np.unique( np.append( all_genes, anno.index.values ) )
    all_genes = np.unique( np.append( all_genes, [string_net.protein1.values, string_net.protein2.values] ) )
    all_genes = np.unique( np.append( all_genes, [op_table.SysName1.values, op_table.SysName2.values] ) )
    all_genes = all_genes[ all_genes != nan ]
    print all_genes.shape

    # gene_regex = get_regex(collect(keys(all_genes)))
    # println(gene_regex)
    
    # ##all_rows::Vector{Int} = [ all_genes[ g ] for g in rownames(x) ] ##pre-define this (not used anymore?)
    
    # ## Get all upstream seqs and their bgFreqs
    # ## TODO: Need to add 0's for k-mers that are NOT in the sequences.
    # ##   Use generate_all_kmers() for that.
    # ## TODO: Need to include vague IUPAC symbols better
    # all_seqs = get_sequences(anno["sysName"].data,anno,genome_seqs,true,op_table,distance_search,false); 
    # all_seqs = all_seqs[ find(all_seqs[:,1].!=""), : ] 
    # all_seqs = filter_sequences( all_seqs, distance_search )

    # all_seqs_scan = get_sequences(anno["sysName"].data,anno,genome_seqs,true,op_table,distance_scan,false); 
    # all_seqs_scan = all_seqs_scan[ find(all_seqs_scan[:,1].!=""), : ]
    # all_seqs_scan = filter_sequences( all_seqs_scan, distance_scan )

    # allSeqs_fname = "./$(organism)/allSeqs.fst"
    # writeFasta( all_seqs_scan, allSeqs_fname ) ## NOTE all_seqs_scan are not used from here on
    
    # ##bgCounts = getBgCounts( [genome_seqs,revComp(genome_seqs)], [0:5], true );
    # ##bgFreqs = getBgFreqs( bgCounts ); ## TODO: bgFreqs are currently not used in MEME-ing OR MAST-ing.

    # all_bgCounts = getBgCounts( all_seqs["seq"].data );
    # all_bgFreqs = getBgFreqs( all_bgCounts );  ## TODO: bgFreqs are currently not used in MEME-ing OR MAST-ing.

    # save_jld( "./$organism/data.jldz", (organism, k_clust, ratios, genome_seqs, anno, op_table, string_net, 
    #                                       allSeqs_fname, all_bgFreqs, all_genes) ) ##, all_rows
    
    # (ratios, genome_seqs, anno, op_table, string_net, ##all_seqs, all_seqs_scan, 
    #  allSeqs_fname, all_bgFreqs, all_genes) ##, all_rows) ##all_bgCounts, 
    return (ratios, genome_seqs, anno, string_net, op_table, all_genes)

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

    good_rows = (x==NA).apply(sum, axis=1) < x.shape[1]/2
    x = x[good_rows]
    print x.shape

    ##for i=1:size(x,1) x[i,:] = stdize_vector(x[i,:]); end ##(x[i,:]-nanmean(x[i,:]))/nansd(x[i,:]); end
    f = lambda x: (x-np.mean(x.dropna())) / np.std(x.dropna())
    x = x.apply(f, axis=1)
    print x.ix[:, :6].head()
    return x

def load_genome(organism):
    ## see http://biopython.org/wiki/SeqIO#Sequence_Input
    ## TODO: allow for multiple genome seqs (e.g. Halo, Sce)
    org_files = np.array( os.listdir('./' + organism + '/') )
    genome_file = './' + organism + '/' + org_files[ np.array( [f.startswith('genome.') for f in org_files] ) ][ 0 ]
    print genome_file

    ##handle = gzip.open('./' + organism + '/' + genome_file)
    ##genome_seqs = SeqIO.to_dict(SeqIO.parse(handle, "fasta")) ## allows >1 sequence too but the files only have one seq so
    ##handle.close()
    genome_seqs = SeqIO.read(gzip.open(genome_file), "fasta") ## this is a python dict
    ## get the biosequence via: genome_seqs.values()[0].seq

    print genome_seqs
    print len(genome_seqs)
    return genome_seqs

def load_annos(organism):
## Load the gene annotations
    org_files = np.array( os.listdir('./' + organism + '/') )
    genomeInfo_file = './' + organism + '/' + org_files[ np.array( [f.startswith('genomeInfo.') for f in org_files] ) ][ 0 ]

    print genomeInfo_file
    ## note x.ix[:,:5].head() prints first 5 cols
    x = pd.read_table(genomeInfo_file, compression='gzip', index_col='sysName') 
    return x

def load_string_net(organism):
## Load the string network
    org_files = np.array( os.listdir('./' + organism + '/') )
    string_file = './' + organism + '/' + org_files[ np.array( [f.startswith('STRING') for f in org_files] ) ][ 0 ]

    print string_file
    string_net = pd.read_table(string_file, compression='gzip', names=['protein1','protein2','weight']) 

    ## Symmetrize it? -- no, seems to already be done.
#   tmp = pd.concat([string_net.ix[:,1], string_net.ix[:,0], string_net.ix[:,2]], axis=1)
#   tmp.columns = tmp.columns[[1,0,2]] ## reorder the columns to the same as string_net
#   string_net = pd.concat( [string_net, tmp] )

    print string_net.shape
    return string_net

def load_op_table(organism):
## Now, require the operons table!
## DONE? (need to verify): Allow for no operons table (e.g. yeast)
    org_files = np.array( os.listdir('./' + organism + '/') )
    op_table = pd.DataFrame()
    try:
        op_file = './' + organism + '/' + org_files[ np.array( [f.startswith('microbesonline_operons_') 
                                                                for f in org_files] ) ][ 0 ]
        print op_file
        op_table = pd.read_table(op_file, compression='gzip') 
        op_table = op_table[ [ "SysName1", "SysName2", "bOp", "pOp" ] ]
        print op_table.shape
    except:
        print "No operons file"
        
    return op_table

## Load the junkey code as text so it can be stored in the run results
# function load_junkey_code(path)
#     files = system(`ls $path/`)
#     files = files[ [ endswith( f, ".jl" ) for f in files ] ]
#     files = files[ files .!= "nanopond-1.9C.jl" ]
#     code = { f => readall( "./$f" ) for f in files }
#     code
# end

# function init_biclusters( x, k_clust, method="kmeans" ) 
#     ## Init via kmeans -- TODO: other methods (including random)
#     ## DONE: use k_means from Clustering package
#     clusters = Dict{Int64,bicluster}()

#     if method == "kmeans" || method == "kmeans+random"
#         xx = x.x;
#         is_nan = isnan(xx)
#         xx[is_nan] = rand(sum(is_nan))*0.1 - 0.05; ## randomize a bit; note this works on the original (global) x!

#         ##km1=kmeansclust( xx, 50, 20 ); ## original
#         ##km1=kmeans( xx', k_clust, 20 ); ## my "optimized" version that pre-removes NAs
#         ##assignments = km1[3]
#         km1=Clustering.kmeans(convert(Matrix{Float64},xx'), k_clust)
#         assignments = km1.assignments
        
#         #     if nprocs() > 1 km2=kmeans2( x.x', k_clust, 20 ); ## parallelized version -- now seems to be broken
#         #     else km1=kmeans( x.x', k_clust, 20 ); ## my "optimized" version that pre-removes NAs
#         #     end
        
#         ## seed each bicluster with rows=output from kmeans and cols=random (1/2 of all cols)
#         for k=1:k_clust
#             rows = findin(assignments, k)
#             if method == "kmeans" 
#                 clusters[k] = bicluster( k, rows, x )
#             elseif method == "kmeans+random" ##         ## add some random rows to each bicluster
#                 clusters[k] = bicluster( k, unique([ rows, randperm(size(xx,1))[1:length(rows)] ]), x )
#             end
#         end
#         xx[is_nan] = NA ## reset the values to NA
#     elseif method == "random"
#         for k=1:k_clust
#             clusters[k] = bicluster( k, unique( randperm(size(x.x,1))[1:20] ), x )
#         end
#     end
#     clusters
# end

# #get_regex( strings::Vector{ASCIIString} ) = get_regex( strings, 2 )

def get_regex( strings, min_ignore=2 ):
    nchar = np.array( [ len(i) for i in strings ] )
    out = ''
    for i in arange( nchar.max() ):
#         d = Dict{Char,Int64}()
#         for j=1:length(strings)
#             if i > nchar[j] continue; end
#             c::Char = strings[j][i]
# 	    d[c] = get(d,c,0) + 1
#         end
#         if isempty(d) break; end
#         tmp = ""
#         for c::Char=collect(keys(d)) if get(d,c,0) > min_ignore tmp = "$tmp$c"; end; end ## ignore values with <=2 occurences
#         if length(tmp) > 1 tmp = "[$tmp]"; end
#         if sum(collect(values(d))) < length(strings) * 0.95 tmp = "$tmp?"; end
#         if tmp != "" out = "$out$tmp"; end
#     end
#     out
# end
