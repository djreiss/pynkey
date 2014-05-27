## DONE: allow parallelization by adding seqs to the bicluster struct, sending biclusters (and allSeqs_fname) 
##    through to all child nodes. Need to make a function that takes a single argument (Dict?) and use pmap().
## DONE: Add n_motifs to the schedule so it ups to 2 after n_iters*3/4 iterations -- NOTE needs to get passed to
##    child processes
## DONE: add min/max_width and search/scan_distance paramters as globals (e.g. to change for yeast)
## DONE: add filter_sequences function to run dust and remove ATGs

##re_meme_all_biclusters_parallel( clusters ) = re_meme_all_biclusters_parallel( clusters, false )

function get_n_motifs()
    global iter, n_iters
    get_n_motifs( iter, n_iters )
end

function get_n_motifs( iter, n_iters )
    n_motifs = ( iter < n_iters * 3 / 4 ) ? 1 : 2 ## n_motifs is 1 in early iterations; 2 in later.
    n_motifs
end

## DONE: only send (and receive) the meme_out and mast_out data (not entire bicluster!) -- a bit faster..?
#re_meme_all_biclusters_parallel( clusters::Dict{Int64,bicluster}, force::Bool ) =
#    re_meme_all_biclusters_parallel( clusters, force, false )

function re_meme_all_biclusters_parallel( clusters::Dict{Int64,bicluster}, force::Bool=false, verbose::Bool=false )
    global k_clust, allSeqs_fname, iter, n_iters, distance_search
    data::Array{Any,1} = []
    n_motifs = get_n_motifs( iter, n_iters )
    r_rownames = rownames(ratios)
    for k in 1:k_clust
        b = clusters[k]
        if ! force && ! b.changed[1] continue; end ## If rows not changed, skip it
        seqs = get_sequences( r_rownames[b.rows] ) 
        seqs = filter_sequences( seqs, distance_search )
        dat = { "k" => b.k, "seqs" => seqs, "allSeqs_fname" => allSeqs_fname, "n_motifs" => n_motifs, "verbose" => verbose };
        push!( data, dat )
    end
    memeOut = pmap( re_meme_bicluster, data ) ## returns an Array{Any,1}
    for i in 1:length(memeOut)
        if memeOut[i][2] == nothing continue; end
        k = memeOut[i][1]
        clusters[k].meme_out = memeOut[i][2]
        clusters[k].mast_out = memeOut[i][3]
    end
    clusters
end

## for parallelization from re_meme_all_clusters_parallel() :
function re_meme_bicluster( x::Dict{Any,Any} ) 
    if x["seqs"] == nothing 
        k = x["k"] ## shouldn't happen!
        warn( "NO GENES! $k" )        
        return( (k, nothing, nothing) )
    end 
    re_meme_bicluster(x["k"], x["seqs"], x["n_motifs"], x["allSeqs_fname"], x["verbose"])
end

function re_meme_all_biclusters( clusters::Dict{Int64,bicluster}, force::Bool=false )
    global k_clust, iter, n_iters
    n_motifs = get_n_motifs( iter, n_iters )
    for k in 1:k_clust
        b = clusters[k]
        if ! force && ! b.changed[1] continue; end ## If rows not changed, skip it
        b = re_meme_bicluster( b, n_motifs )
        clusters[k] = b
    end
    clusters
end

function re_meme_bicluster( b::bicluster, n_motifs::Int=2 )
    global ratios, allSeqs_fname, distance_search
    seqs = get_sequences( rownames(ratios)[b.rows] ) 
    seqs = filter_sequences( seqs, distance_search )
    (k, meme_out, mast_out) = re_meme_bicluster( b.k, seqs, n_motifs, allSeqs_fname, false )
    b.meme_out = meme_out
    b.mast_out = mast_out
    b
end

function re_meme_bicluster( k::Int64, seqs::DataFrame, ##Matrix{ASCIIString}, 
                           n_motifs::Int, allSeqs_fname::ASCIIString, verbose::Bool=false )
    #@time gibbs_out = gibbs_site_sampler(seqs["seq"].data)     ## run gibbs sampler on most "flagellar-enriched" cluster
    #@time gibbs_out2 = gibbs_site_sampler(seqs, gibbs_out["pssm"])

    meme_out::Vector{ASCIIString} = []
    mast_out::DataFrame = DataFrame()
    try
        ns = size(seqs,1)
        println("IN HERE: MOT $n_motifs $k $ns")
        meme_out = do_meme(seqs,n_motifs,verbose)
        if meme_out != "" # NOTE for now if meme fails (most often b/c too few sequences) we just keep the previous run
            try
                mast_out = do_mast(meme_out,allSeqs_fname,false,verbose)
            catch x
                warn( "ERROR RUNNING MAST FOR BICLUSTER $k" )
                println(x)
            end
        else
            warn( "TOO FEW SEQUENCES TO MEME FOR BICLUSTER $k" )
        end
    catch x
        warn( "TOO FEW SEQUENCES TO MEME FOR BICLUSTER $k" )
        ##println(x)
    end
    (k, meme_out, mast_out)
end

## This is parallelizable because it sends the seqs to be searched to each node
function do_meme(seqs::DataFrame, n_motifs::Int=2, verbose::Bool=false) 
    global iter, n_iter, motif_width_range
    seqs = seqs[ seqs["seq"].data .!= "", : ]
    seqs = seqs[ ! duplicated(seqs["seq"].data), : ]
    if size(seqs,2) <= 2 return( "" ); end

    (seqs_fname,str) = mktemp()
    writeFasta( seqs, seqs_fname )
    
    cmd = `./progs/meme $seqs_fname -time 600 -dna -revcomp -maxsize 9999999 -nmotifs $n_motifs -evt 999999 -minw $(motif_width_range[1]) -maxw $(motif_width_range[2]) -mod zoops -nostatus -text`
    if verbose println( cmd ); end
    
    memeOut = system(cmd)
    rm( seqs_fname )

    #if verbose print(join(memeOut, '\n')); end
    memeOut
end

function do_mast(memeOut, verbose::Bool=true)
    global allSeqs_fname
    do_mast(memeOut, allSeqs_fname, false, verbose)
end
                        
## this IS parallelizable b/c it just uses the allSeqs_fname file for reading, which was already created
function do_mast(memeOut, allSeqs_fname, get_allHitsTab=false, verbose=true) 
    (memeOutFname,str) = mktemp()
    for i=memeOut write( str, "$i\n" ); end
    close( str )

    ## Not used right now but TODO: convert output to DataFrame (currently an Any array)
    if get_allHitsTab
        ## First get the table of sites (multiple hits per sequence) -- using "-hits_list"
        cmd = `./progs/mast $memeOutFname -d $allSeqs_fname -nostatus -stdout -text -brief -ev 999999 -mev 999999 -mt 1.0 -seqp -remcorr -hit_list`
        if verbose println(cmd); end
        mastOut = system(cmd)

        tmpHitsTab = Array(Any, length(mastOut)-4, 5)
        i = 1
        for j=3:length(mastOut)-2
            tmp = split( mastOut[j], r"\s+" )
            tmpHitsTab[ i,1 ] = tmp[1]
            tmpHitsTab[ i,2 ] = parseint(tmp[2])
            tmpHitsTab[ i,3 ] = parseint(tmp[3])
            tmpHitsTab[ i,4 ] = parseint(tmp[4])
            tmpHitsTab[ i,5 ] = parsefloat(tmp[5])
            #tmpHitsTab[ i,6 ] = parsefloat(tmp[6])
            i += 1
        end

        mastAllHitsTab = DataFrame( { "Gene" => convert( Array{ASCIIString}, tmpHitsTab[ :,1 ] ),
                                      "Mot" => convert( Array{Int32}, tmpHitsTab[ :,2 ] ),
                                      "Start" => convert( Array{Int64}, tmpHitsTab[ :,3 ] ),
                                      "End" => convert( Array{Int64}, tmpHitsTab[ :,4 ] ),
                                      "E-value" => convert( Array{Float32}, tmpHitsTab[ :,5 ] ) } )
    end ## get_allHitsTab

## This gets the table of sequence e-values (one value per sequence) -- need to parse it out
## Table starts 2 lines below 
## SEQUENCE NAME                      DESCRIPTION                   E-VALUE  LENGTH
## and last line is 2 lines above 
## ********************************************************************************
    cmd = `./progs/mast $memeOutFname -d $allSeqs_fname -nostatus -stdout -text -brief -ev 999999 -mev 999999 -mt 0.99 -seqp -remcorr` #-b
    if verbose println(cmd); end
    mastOut = system(cmd)

## We need P-values, not E-values, so we turned off the '-b' flag and need to parse the 3rd section (here is the 2nd section)
## We can get the e-vals from the 3rd section too, so let's use that
    # startInd = find( mastOut .== "SEQUENCE NAME                      DESCRIPTION                   E-VALUE  LENGTH" )[ 1 ] + 2
    # endInd = find( mastOut .== "********************************************************************************" )
    # endInd = endInd[ endInd .> startInd ][ 1 ] - 2
    # mo = mastOut[ startInd:endInd ]
    # ##mastOutTab = Array(Any, length(mastOut), 2)
    # mo_seq = Array( ASCIIString, length(mo) )
    # mo_ev = Array( Float32, length(mo) )
    # for i=1:length(mo)
    #     tmp = split( mo[i], r"\s+" )
    #     mo_seq[ i ] = tmp[1]
    #     mo_ev[ i ] = parsefloat(tmp[2])
    # end
    # mastOutTab = DataFrame( { "Gene" => mo_seq, "E-value" => mo_ev } )

## OK, now parse the 3rd section
    reg = r"COMBINED P-VALUE"
    lines = find( [ ismatch( reg, mastOut[ i ] ) for i=1:length(mastOut) ] .== true )
    genes = convert( Array{ASCIIString}, mastOut[ lines-2 ] )
    reg = r"COMBINED P-VALUE = (.+)(\s+E-VALUE)"
    tmp = [ match(r"COMBINED P-VALUE = (.+)(\s+E-VALUE = (.*))", mastOut[i]).captures[[1,3]] for i=lines ]
    pvals = float32( [ parsefloat( strip(tmp[i][1]) ) for i=1:length(tmp) ] )
    evals = float32( [ parsefloat( strip(tmp[i][2]) ) for i=1:length(tmp) ] )
    mastOutTab = DataFrame( { "Gene" => genes, "E-value" => evals, "P-value" => pvals } )

    rm( memeOutFname )
    #rm( allSeqs_fname )
    if get_allHitsTab return(mastOutTab, mastAllHitsTab)
    else return(mastOutTab)
    end
end

