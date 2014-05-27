## TODO: handle fasta files w/ multiple entries. Right now, concatenates them all into one seq.
## TODO: use BioSeq pkg (can use sub- and file-based arrays, and DNA2seq which is 2 bits per NT!) -- see bottom for sample code 
##    NOTE that DNA2seq can't handle IUPAC symbols but file-based can! (and is just as fast!)
## TODO: Need to add remove-ATG step, and also running dust on sequences.

#get_sequences( genes::Vector{ASCIIString}, op_shift::Bool, distance::Vector{Int64} ) =
#                           get_sequences( genes, op_shift, distance, false )

#get_sequences( genes::Vector{ASCIIString} ) = get_sequences( genes, true, distance_search )

## This version of the function only works after anno/genome_seq/op_table have been set to a global var.
function get_sequences( genes::Vector{ASCIIString}, op_shift::Bool=true, distance::Vector{Int64}=distance_search, 
debug::Bool=false )
    global anno, genome_seqs, op_table;
    get_sequences( genes, anno, genome_seqs, op_shift, op_table, distance, debug )
end

function get_sequences( genes::Vector{ASCIIString}, anno::DataFrame, genome_seqs::DataFrame, 
                       op_shift::Bool, op_table::DataFrame, distance::Vector{Int64}, debug::Bool=false )
    ##genes = genes[ findin( genes, rownames(anno) ) ]
    seqs::Vector{DataFrame} = Array(DataFrame, length(genes)) ##Matrix{ASCIIString} = Array(ASCIIString,length(genes),3)
    ind = 0
    for gene=genes
        #println(gene," ",distance)
        anno_ind = findfirst( anno["sysName"].data, gene )
        upstream_gene = gene
        if anno_ind == 0 ##! any(anno["sysName"].data .== gene) 
            ##seqs[ ind += 1,: ] = [gene, upstream_gene, ""]
            seqs[ ind += 1 ] = DataFrame( {"gene"=>[gene], "upstream_gene"=>[upstream_gene], "seq"=>[""]} )
            continue
        end
        if op_shift && size(op_table,1) > 0
            strand = anno["strand"][ anno_ind ] ##findin(anno["sysName"].data, gene)[1] ]  ##anno["sysName"].data.==gene,"strand"][1]
            upLabel = strand == "+" ? "SysName1" : "SysName2"
            currLabel = strand == "+" ? "SysName2" : "SysName1"
            row = findfirst( op_table[ currLabel ], gene ) ##find( op_table[ currLabel ] .== gene )
            while row > 0 && op_table[ row, "pOp" ] >= 0.8
                upstream_gene = op_table[ row, upLabel ]
                ##println(gene," ",row," ",upstream_gene," ",op_table[ row[1], "pOp" ])
                row = findfirst( op_table[ currLabel ], upstream_gene ) ##find( op_table[ currLabel ] .== upstream_gene )
            end 
        end
        tmp_ind = findfirst( anno["sysName"].data, upstream_gene )
        strnd = anno[ "strand" ][ tmp_ind ] ##anno[ anno["sysName"].data.==upstream_gene, "strand" ][1]
        strt = anno[ "start" ][ tmp_ind ] ##anno[ anno["sysName"].data.==upstream_gene, "start" ][1]
        rng = strnd == "+" ? strt + distance : strt - reverse(distance)
        genome_seq = get_genome_seq( genome_seqs, anno[ "scaffoldId" ][ tmp_ind ] )
        if rng[2] > length(genome_seq) rng[2] = length(genome_seq); end
        seq = genome_seq[ rng[1]:rng[2] ]
        if debug println( "$gene $upstream_gene $strnd $strt $rng" ); end
        ##println("$gene $upstream_gene $strnd $strt $rng $seq")
        if strnd == "-" seq = revComp(seq); end
        ##seqs[ ind+=1,: ] = [gene, upstream_gene, seq]
        seqs[ ind+=1 ] = DataFrame( {"gene"=>[gene], "upstream_gene"=>[upstream_gene], "seq"=>[seq]} )
    end
    ##seqs
    out = rbind(seqs)
    out
end

function filter_sequences( seqs::DataFrame, distance=distance_search, remove_repeats=true, remove_atgs=true )
    seqs = seqs[ seqs["seq"] .!= "", : ] ## remove all empty sequences; they cause heartburn.

    if remove_repeats == true ##&& length( grep( "NNNNNN", seqs ) ) <= 1
       ##if verbose println( "Removing low-complexity regions from sequences.\n" ); end
       fname = tempname()
       writeFasta( seqs, fname )
       tmp = system( `./progs/dust $fname` ) ## note dust as an optional numeric parameter -- what does it mean?
       rm( fname )
       (fname,io) = mktemp()
       for line in tmp write(io, "$(line)\n"); end
       close( io )
       seqs_new = readFastaDNA( fname )
       rm( fname )
       if size(seqs,1) != size(seqs_new,1) 
           warn( "Remove low complexity failed - skipping!" )
       elseif all(seqs["gene"] .== seqs_new["gene"])
           seqs_new["upstream_gene"] = seqs["upstream_gene"]
           seqs = seqs_new
       end
    end

    if remove_atgs && any( distance .< 0 ) 
        nnnn = convert(Vector{Uint8}, "NNNN")
        for i in 1:size(seqs,1)
            ss = convert(Vector{Uint8}, seqs["seq"][i])
            ss[ (distance[2]+1):(distance[2]+4) ] = nnnn  ## Mask out ATGs (why 4 instead of 3?)
            seqs["seq"][i] = convert(ASCIIString, ss)
        end
    end

    seqs
end

function get_genome_seq( genome_seqs::DataFrame, scaffoldId::Int64 )
    sid = "$scaffoldId"
    ind = find( [beginswith(genome_seqs["gene"].data[i], sid) for i=1:size(genome_seqs,1)] )[ 1 ]
    genome_seqs[ "seq" ].data[ ind ]
end

## Create DataFrame for multiple sequences; same format as get_sequences(bicluster)
function readFastaDNA( fname ) 
    str = open( fname )
    seq::ASCIIString = readall( str );
    close( str );
    seqs = split( seq, '>' );
    first_return = [search(seqs[i], '\n') for i=2:length(seqs)] ## first seq is ""
    seq_names = convert( Vector{ASCIIString}, [seqs[i][1:(first_return[i-1]-1)] for i=2:length(seqs)] ); 
    seqs = [seqs[i][(first_return[i-1]+1):end] for i=2:length(seqs)];
    seqs = convert( Vector{ASCIIString}, [uppercase( replace( replace( seqs[i], r">.*\n", "" ), '\n', "" ) ) for i=1:length(seqs)] );
    out = DataFrame( {"gene"=>seq_names, "upstream_gene"=>fill("", length(seqs)), "seq"=>seqs} );
    out
end

function writeFasta( seqs, fname ) ## Assumes seqs in DataFrame format of get_sequences()
    str = open( fname, "w" )
    for i=1:size(seqs,1)
        gene = seqs["gene"].data[i]
        seq = seqs["seq"].data[i]
        write( str, ">$gene\n" )
        write( str, "$seq\n" )
    end
    close( str )
end

function revComp( seq::ASCIIString )
    out = Array( Uint8, length( seq ) ) ## uninitialized is faster?
    j = length( seq ) + 1
    seq = uppercase( seq )
    for i = 1:length( seq )
        c = seq[ i ]
        out[ j -= 1 ] = c == 'G' ? 'C' : ( c == 'C' ? 'G' : ( c == 'T' ? 'A' : ( c == 'A' ? 'T' : c ) ) )
    end
    ASCIIString( out );
end

const DNA_letters = [ 'G', 'A', 'T', 'C' ];
const DNA_letter_lookup = {'G'=>1, 'A'=>2, 'T'=>3, 'C'=>4 };

function read_iupac()
    fname = "./IUPAC-dna.txt"
    str = open( fname )
    lines = split( readall( str ), '\n' )
    close( str )
    lines = [ replace(l," ","") for l=lines]
    lines = [ replace(l,r"\(.*\)","") for l=lines ]
    lines = lines[ lines .!= "" ]
    d = Dict()
    for l in lines
        if any(DNA_letters .== l[1]) d[l[1]] = l[1]; 
        elseif l[1] == 'U' d[l[1]] = l[1]; 
        else d[l[1]] = l[3:end]; end
    end
    d
end

const iupac_dict = read_iupac();

# function collect_add_dicts( a, b ) ## Add up all elements into output dict
#     for (k,v) in b a[k] = v + b[k]; end
#     a
# end

#@profile begin
## Note this will not have an entry for subseqs that don't exist in the sequence
## Only does it for a SINGLE "order" value and single strand - see the default getBgCounts() to see standard usage.
## Use in conjunction with generate_all_kmers in order to get seqs w/ count=0
## Taken from my seqTest.jl
function getBgCounts( seqs::Array{ASCIIString,1}, order::Array{Int64,1}=[0:5], verbose::Bool=false ) 
    ss::ASCIIString = "";
    d = Dict{ASCIIString,Int64}() ## Specifying the types speeds it up about 35%
    for seq=[seqs,[revComp(s) for s=seqs]] ##seqs
        seq = uppercase( seq )
        for j=order
            if verbose println(j); end
            for i=1:length(seq)-j ## for loop is about 5% faster than while loop
	        ss = seq[i:(i+j)]
	        d[ss] = get(d,ss,0) + 1
            end
        end
    end
    d
end

##getBgCounts( seqs::Array{ASCIIString,1}, order::Int64 ) = 
##             getBgCounts( [seqs,[revComp(s) for s=seqs]], [order] );
# getBgCounts( seqs::Array{ASCIIString,1}, order::Array{Int64,1} ) = 
#                getBgCounts( [seqs,[revComp(s) for s=seqs]], order ); ##, false );
# getBgCounts( seqs::Array{ASCIIString,1} ) = 
#                getBgCounts( [seqs,[revComp(s) for s=seqs]] );
# getBgCounts( seq::ASCIIString ) = getBgCounts( [seq] );

#bgCounts = getBgCounts( genome_seqs );
#end # profile

## Convert counts dictionary into frequencies; divide counts for each k-mer by the total for that given k.
function getBgFreqs( bgCounts::Dict{ASCIIString,Int64} )
    k = collect(keys( bgCounts ))
    nc = [ length(k[i]) for i=1:length(k) ]
    d = Dict{ASCIIString,Float64}()
    for i=1:maximum(nc)
        ks = k[find(nc.==i)]
        tot = sum( [bgCounts[j] for j=ks] )
        for j=ks d[j] = bgCounts[j] / tot; end
    end
    d
end

## Couldn't do this in R! Generate all sequences with given length.
## Allow for up to n_ns N's as well.
function generate_all_kmers( len::Int64, n_ns::Int64=0, letters::Array{Char,1}=DNA_letters )
    function append_nucs( d::Array{ASCIIString}, letters::Array{Char,1} )
        out = Array(ASCIIString, length(d)*length(letters))
        ind::Int64 = 0
        for i=d for j=letters out[ind+=1] = strcat(i,j); end; end ## note can also use i*j
        out
    end

    if n_ns > 0 letters = [ letters, 'N' ]; end ## Add N to letters; just need 1 copy.
    d = Array(ASCIIString, length(letters))
    i=0; for j=letters d[i+=1] = "$j"; end
    for i = 1:len-1 d = append_nucs( d, letters ); end
    if n_ns > 0   ## Get rid of k-mers with > n_ns Ns.
        n_count = [ sum(chars(d[i]) .== 'N') for i=1:length(d) ]
        d = d[ find( n_count .<= n_ns ) ]
    end
        
    d
end

#generate_all_kmers( len::Int64 ) = generate_all_kmers( len, 0, DNA_letters )
#generate_all_kmers( len::Int64, n_ns::Int64 ) = generate_all_kmers( len, n_ns, DNA_letters )

if false
    seq = "GATCATGCATGTATGCTACGTGCGCGGGTACGTATATGATGCTATTATCGTAGCTACGTAGCTAGCTAGCTACAGTCGATCGATTGAC"
    seq2bit = dna2seq(seq)
    d=Dict{BitArray{1},Int64}()
    d[[seq2bit[1:5].b1,seq2bit[1:5].b2]] = 1 ## can count up k-mers same as with ASCIIString

    ## file-based - can handle IUPAC codes
    seqnt = nt(seq)
    Astream=open("test.mmap", "a+")
    A=mmap_array(Nucleotide, size(seqnt), Atream)
    A[1:length(seqnt)] = seqnt
    close(Astream)

    ## try this! file-based, don't even have to read it in!
    len = filesize("Hpy/genome.85962.txt")
    seqStream = open("Hpy/genome.85962.txt", "r")
    seq = mmap_array(Nucleotide, (len,), seqStream)
    println(convert(ASCIIString,seq[500:505]))
end