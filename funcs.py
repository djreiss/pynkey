function run_junkey()
    global organism, iter, n_iters, clusters, stats_df

    n_no_improvements = 0
    for i=iter:n_iters
        iter = i
        (clusters, n_improvements, stats_tmp) = do_floc( clusters );
        println( @sprintf( "%.3f", (time() - startTime)/60 ), " minutes since initialization" )
        stats_df = rbind( stats_df, stats_tmp )
        write_table( "output/$(organism)_stats.tsv", stats_df )
        if isfile( "DO_SAVE" )  ## save cluster info for temporary examination of clusters (via Rscripts/clusters.R)
            warn( "Writing out clusters to output/$(organism)_clusters.tsv" )
            clusters_tab = clusters_to_dataFrame(clusters);
            write_table("output/$(organism)_clusters.tsv", clusters_tab)
        end
        if n_improvements <= 0 n_no_improvements += 1 else n_no_improvements = 0; end
        if iter > n_iters/2 && n_no_improvements > 5 break; end
    end
end

function copy_cluster( cc::bicluster, full=true, everything=true )
    out = cc
    if everything ## make copies of everything
        out = bicluster( cc.k, copy(cc.rows), copy(cc.cols), cc.var, cc.resid, cc.dens_string, cc.meanp_meme,
                        copy(cc.scores_r), copy(cc.scores_c), copy(cc.scores_n), copy(cc.scores_m), 
                        copy(cc.meme_out), copy(cc.mast_out), cc.changed )
    elseif full ## Make copies of rows/cols so if they get updated in the original, they won't change here.
        out = bicluster( cc.k, copy(cc.rows), copy(cc.cols), cc.var, cc.resid, cc.dens_string, cc.meanp_meme, 
                        cc.scores_r, cc.scores_c, cc.scores_n, cc.scores_m, 
                        cc.meme_out, cc.mast_out, cc.changed )
    else ## This is a shallow copy - changes to rows/cols in this cluster will mirror those in cc
        out = bicluster( cc.k, cc.rows, cc.cols, cc.var, cc.resid, cc.dens_string, cc.meanp_meme, 
                        cc.scores_r, cc.scores_c, cc.scores_n, cc.scores_m, 
                        cc.meme_out, cc.mast_out, cc.changed )
                   
    end
    out
end

function copy_clusters( clusters::Dict{Int64,bicluster}, full=true, everything=true )
    global k_clust
    new_clusters = Dict{Int64,bicluster}() ## make a copy for updating
    for k=1:k_clust new_clusters[k] = copy_cluster( clusters[k], full, everything ); end
    new_clusters
end

function bicluster_residue{T}( b::bicluster, x::NamedMatrix{T} )
    rats = x.x[ b.rows, b.cols ]
    if ndims( rats ) < 2 || size( rats, 1 ) <= 1 || size( rats, 2 ) <= 1 || mean( isnan( rats ) ) > 0.95 
        warn( "COULD NOT COMPUTE RESIDUE FOR BICLUSTER $(b.k)" )
        return( convert(T,1.0) )
    end

    d_rows = rowmeans( rats )
    d_cols = colmeans( rats )
    d_all = nanmean( d_rows )
    rats += d_all - outer( d_rows, d_cols, + )

    average_r = nanmean( abs( rats ) )
    average_r
end

function bicluster_var{T}( b::bicluster, x::NamedMatrix{T} )
    const var_add = 0.1
    x1::Matrix{T}=x.x[b.rows,b.cols] ## get the bicluster's submatrix of the data
    mn::Vector{T}=colmeans(x1) ## subtract bicluster mean profile
    for i=1:size(x1,2) x1[:,i] -= mn[i]; end
    nanvar(x1) / (nanvar(mn) + var_add)
end

bicluster_volume( b::bicluster ) = float32( length(b.rows) * length(b.cols) )

function re_seed_bicluster_if_necessary( clust::bicluster )
    const min_rows = 3; const max_rows = 80
    if length(clust.rows) < min_rows ||  ## Need to fill in "empty" clusters, mostly because var is non-defined, Also because meme-ing only works with >2 sequences
        length(clust.rows) > max_rows   ### Squash clusters that get way too big (DONE: just remove some to bring it down to max_rows)
        nr = length(clust.rows)
        warn( "RESEEDING BICLUSTER $(clust.k) ($nr)" )
    ## DONE: add rows that are preferentially in few (or no) other clusters -- sample from get_cluster_row_counts(clusters)
    ##clust.rows = unique([clust.rows, [Distributions.sample([1:nrow(ratios.x)]) for i=1:5]]) ## add 5 random rows
        counts_g::Vector{Float32} = float32(get_cluster_row_counts(clusters))
        if length(clust.rows) < min_rows 
            counts_g = maximum(counts_g) + 0.01 - counts_g
            counts_g = counts_g / maximum(counts_g)
            clust.rows = unique( my_sampler( [1:length(counts_g)], 5, counts_g ) )
            #clust.rows = unique([clust.rows, [Distributions.sample([1:length(counts_g)], counts_g) for i=1:5]]) 
            clust.cols = unique([clust.cols, [Distributions.sample([1:ncol(ratios.x)]) for 
                                              i=1:(div(ncol(ratios.x),3)-length(clust.cols))]])  ## add x/3 random cols
        elseif length(clust.rows) > max_rows
            ##clust.rows = unique([Distributions.sample([1:length(counts_g)], counts_g) for i=1:max_rows]) ## remove some genes, include mostly from genes in few clusters
            ##clust.rows = unique( my_sampler( [1:length(counts_g)], max_rows, counts_g ) ) ## remove some genes, include mostly from genes in few clusters
            counts_g = ( counts_g + 0.01 ) / ( maximum(counts_g) + 0.01 )
            tmp_rows = unique( my_sampler( clust.rows, length(clust.rows)-max_rows+1, counts_g[clust.rows] ) ) ## remove some genes, mostly that are in lots of clusters
            ##println("REMOVING: ",tmp_rows," ",get_cluster_row_counts(clusters)[tmp_rows])
            clust.rows = clust.rows[ ! in(clust.rows, tmp_rows) ]
        end
    end
    clust
end

function re_seed_all_clusters_if_necessary( clusters::Dict{Int64,bicluster} )
    global k_clust
    for k in 1:k_clust clusters[k] = re_seed_bicluster_if_necessary(clusters[k]); end
    clusters
end

#fill_cluster_scores(clust::bicluster) = fill_cluster_scores(clust, false)

function fill_cluster_scores( clust::bicluster, force=false ) 
    global iter
    fill_cluster_scores(clust, iter, force, false)
end

function fill_cluster_scores( clust::bicluster, iter::Int64, force=false, verbose=false )
    global ratios, string_net
    ##clust = re_seed_bicluster_if_necessary(clust)
    if verbose println(clust.k, " ", length(clust.rows), " ", length(clust.cols)); end
    (weight_r, weight_n, weight_m, weight_c, weight_v) = get_score_weights(iter)
    #println("HERE5 $weight_r $(clust.changed)")
    if force || ( ( weight_r + weight_c > 0 ) && ( sum(clust.changed) > 0 ) )  
        try
            clust = get_cluster_expr_rowcol_scores( clust, ratios ) ## Actually need to do it for rows too, even if only cols changed
        catch x
            warn( "ERROR WITH ROWCOL SCORES!" )
            println( x )
        end
    end
    #println("HERE6")
    if force || clust.changed[1] ## rows only
        if abs(weight_n) > 0 
            try
                clust = get_cluster_network_row_scores( clust, string_net )
            catch x
                warn( "ERROR WITH NETWORK SCORES!" )
                println( x )
            end
        end
        #println("HERE6a")
        if weight_m > 0 
            try
                clust = get_cluster_meme_row_scores( clust )
            catch x
                warn( "ERROR WITH MEME SCORES!" )
                println( x )
            end
        end
    end 
    #println("HERE7")
    clust.changed[1] = clust.changed[2] = false
    #println("HERE8")
    clust
end

#fill_all_cluster_scores(clusters::Dict{Int64,bicluster}) = fill_all_cluster_scores(clusters, false, false)

## DONE: if anything needs to be parallelized, it's this function. (At least the meme-ing.) -- see re_meme_all_biclusters_parallel
## DONE: Also use the bicluster.changed field consistently and don't re-calc for biclusters that haven't changed.
function fill_all_cluster_scores( clusters::Dict{Int64,bicluster}, force=false, verbose=false )
    global k_clust
    clusters = re_seed_all_clusters_if_necessary(clusters)    
    ## Get all the scores for adding/removing rows/cols from each bicluster
    for k=1:k_clust ## DONE: good possibility for parallelization? Note, we have to send over all the data
        clust = clusters[k]     ## (including ratios, networks) to children.
        if ! force && ! any(clust.changed) continue; end
        if verbose println(k, " ", length(clust.rows), " ", length(clust.cols)); end
        clust = fill_cluster_scores(clust, force)
        clusters[ k ] = clust
    end
    clusters
end

function pre_load_child_nodes() 
    global organism, k_clust, ratios, string_net, all_genes, n_iters
    ## Need to pre-send all data to children -- I dont know how to send a variable in the @everywhere call so use a fixed filename
    ## This only needs to be done once, doing it multiple times causes a big memory leak on the children.
    println( "Sending data to child nodes" )
    fname = "tmp_send_data.jldz"
    save_jld( fname, (organism, k_clust, ratios, string_net, all_genes, n_iters, motif_width_range, distance_search, distance_scan) ) 
                      #genome_seq, anno, op_table, allSeqs_fname, all_bgFreqs, all_genes, all_rows, iter, clusters) )
    @everywhere (organism, k_clust, ratios, string_net, all_genes, n_iters, motif_width_range, distance_search, distance_scan) = load_jld("tmp_send_data.jldz")
    sleep( 5 ) ## make sure children have time to load the data before the file is removed ???
    rm( fname )
end

# function my_pmap(f, lsts...)
#     np = nprocs()
#     n = length(lsts[1])
#     results = cell(n)
#     i = 1
#     # function to produce the next work item from the queue. 
#     # in this case it's just an index. 
#     next_idx() = (idx=i; i+=1; idx)
#     @sync begin
#         for p=1:np
#             if p != myid() || np == 1
#                 @spawnat myid() begin
#                     while true
#                         idx = next_idx()
#                         if idx > n
#                             break
#                         end
#                         #println("HEREXXX: $idx")
#                         results[idx] = remotecall_fetch(p, f,
#                                                         map(L->L[idx], lsts)...)
#                     end
#                 end
#             end
#         end
#     end
#     results
# end

#fill_all_cluster_scores_parallel(clusters::Dict{Int64,bicluster}) = fill_all_cluster_scores_parallel(clusters, false, false)

function fill_all_cluster_scores_parallel( clusters::Dict{Int64,bicluster}, force=false, verbose=false )
    global k_clust, iter
    clusters = re_seed_all_clusters_if_necessary(clusters) ## do this first since it relies on a global "clusters"
    data::Array{Any,1} = []
    for k in 1:k_clust
        if ! force && ! any(clusters[k].changed) continue; end ## If bicluster not changed, skip it
        b = bicluster( k, clusters[k].rows, clusters[k].cols ) ## copy the cluster but only send over the necessary parts
        b.mast_out = clusters[k].mast_out
        if ! force b.changed = clusters[k].changed
        else b.changed[1] = b.changed[2] = true; end
        dat = { "iter" => iter, "k" => b.k, "biclust" => b, "verbose" => verbose, "force" => force }
        push!( data, dat )
    end
    new_clusts = pmap( fill_cluster_scores, data ) ## returns an Array{Any,1}
    @everywhere gc() ## force a gc -- does it help?
    for i in 1:length(new_clusts)
        b = new_clusts[i]  ## new (updated) cluster
        k = b.k
        bb = clusters[k] ## original cluster
        b.meme_out = bb.meme_out ## this did not get updated in b so just copy it over
        if ! force ## THIS STUFF BELOW IS BECAUSE WE CALL fill_cluster_scores(x::Dict) with force=false
            if ! any(bb.changed) ## rows AND cols not changed so scores_r, scores_c were NOT updated in b
                b.scores_r = bb.scores_r
                b.scores_c = bb.scores_c
                b.resid = bb.resid; b.var = bb.var
            end
            if ! bb.changed[1] ## rows not changed so scores_r, scores_n were NOT updated in b
                b.scores_m = bb.scores_m; b.meanp_meme = bb.meanp_meme
                b.scores_n = bb.scores_n; b.dens_string = bb.dens_string
            end
        end
        clusters[k] = b 
    end
    clusters
end

function fill_cluster_scores( x::Dict{Any,Any} ) 
    k = x["k"]
    iter = x["iter"]
    b = x["biclust"]
    force = x["force"]
    verbose = x["verbose"]
    out = fill_cluster_scores(b, iter, force, verbose)   ### note changed force to false; this may be buggy
    out
end

## number of clusters each gene is in - need to compute only once over all clusters
function get_cluster_row_counts( clusters::Dict{Int64,bicluster} )
    global ratios
    counts::Vector{Int32} = []
    if length(clusters[1].scores_r) > 0 counts = zeros(Int32, length(clusters[1].scores_r)) ## beginning - empty clusters
    else counts = zeros(Int32, size(ratios.x,1)); end
    for k in 1:length(clusters) for r in clusters[k].rows counts[r] += 1; end; end
    counts
end

## Use edge density score: sum of edge weights (=number of edges for all weights =1) / number of nodes^2
## This is the faster version that uses SubDataFrames

function bicluster_network_density( b::bicluster, network::DataFrame )
    global ratios;
    ## Note that right now rows in the bicluster is the index into the rownames of the expression matrix!!! TODO: change this!
    rows::Vector{ASCIIString} = rownames(ratios)[b.rows] ##keys( all_genes )[ in( values( all_genes ), rows ) ]
    bicluster_network_density( rows, network )
end

function bicluster_network_density( rows::Vector{ASCIIString}, network::DataFrame )
    net1 = sub( network, findin( network["x1"], rows ) )
    net2 = sub( net1, findin( net1["x2"], rows ) )
    dens = sum( net2["x3"] ) / length(rows)^2 ## Already symmetrized, need to decrease count by 1/2
    log10( dens+1e-9 )
end

function bicluster_meme_pval( b::bicluster )
    global ratios
    if size(b.mast_out,1) <= 0 return( NA ); end

    rows = b.rows
    genes = rownames(ratios)[rows]
    ## TODO: don't only use values for genes that are in the ratios matrix.
    df = sub( b.mast_out, in( b.mast_out["Gene"].data, genes ) )
    mn = nanmean( log10( df["P-value"].data ) )
    mn
end

#clusters_w_func( func, clusters ) = clusters_w_func( func, clusters, 1 )

## Find the flagellar cluster, whew!!!
function clusters_w_func( func::ASCIIString, clusters, n_best=1 )
    global anno, ratios
    inds = findn([ismatch(Regex(func),anno["desc"][i]) for i=1:size(anno,1)] .== true)
    r_rownames = rownames(ratios)
    inds2 = int32([in(r_rownames, anno["sysName"][i]) ? 
                   findn(r_rownames .== anno["sysName"][i])[1] : 0 for i=inds])
    inds2 = inds2[ inds2 .!= 0 ]

    ord = sortperm( int64([length(findin(clusters[k].rows,int64(inds2))) for k=1:length(clusters)]) )
    kInds = ord[ (end-n_best+1):end ] ## The n_best clusters w/ the most genes annotated with "flagell"
    ##kInd = findmax(int64([length(findin(clusters[k].rows,int64(inds2))) for k=1:length(clusters)]))[ 2 ]  
    
    for kInd in kInds
        genes = r_rownames[clusters[kInd].rows] ##rows]
        ##println(genes) ## print the genes
        genes = genes[ findin(genes, anno["sysName"].data) ]
        println(kInd, "\n", anno[in(anno["sysName"].data,genes),["sysName","desc"]])
    end
    kInds
end

function print_cluster_stats( clusters::Dict{Int64,bicluster} )
    global k_clust, iter, startTime
    time_elapsed = (time() - startTime)/60
    ## To get memory used -- DONE: add this to the stats_df; TODO: get this during MEME running (that's when max RAM is used) 
    ## ps -U dreiss --no-headers -o rss -o comm | grep -E 'julia|meme|mast' | awk '{print $1}' | ( tr '\n' + ; echo 0 ) | bc
    if OS_NAME != :Darwin
        tmp = 0
        try ## this command doesnt work on osiris although it works on fossil... weird
            tmp = split( readall(`ps -U dreiss --no-headers -o rss -o comm` | `grep -E 'julia|meme|mast'` | `awk '{print $1}'`), '\n' )
            tmp = sum([ parse_int(tmp[i]) for i=1:(length(tmp)-1) ])
        catch
            tmp = 0
        end
    else
        tmp = split( readall(`ps -U dreiss -o rss -o comm` | `grep -E 'julia|meme|mast'` | `awk '{print $1}'`), '\n' )
        tmp = sum([ parse_int(tmp[i]) for i=1:(length(tmp)-1) ])
    end
    (weight_r, weight_n, weight_m, weight_c, weight_v) = get_score_weights()
    out_df = DataFrame( { "iter" => iter, "time" => time_elapsed, "mem_used" => tmp, 
                         "r0" => weight_r, "n0" => weight_n, "m0" => weight_m, 
                         "c0" => weight_c, "v0" => weight_v } )
    tmp = float32([length(clusters[k].rows) for k=1:k_clust])
    out_df["ROWS"] = nanmean(tmp)
    println( "ROWS: ", out_df["ROWS"][1], " +/- ", nansd(tmp) )
    tmp = float32([length(clusters[k].cols) for k=1:k_clust])
    out_df["COLS"] = nanmean(tmp)
    println( "COLS: ", out_df["COLS"][1], " +/- ", nansd(tmp) )
    tmp = float32([clusters[k].resid for k=1:k_clust])
    out_df["RESID"] = nanmean(tmp)
    println( "MEAN RESID: ", out_df["RESID"][1], " +/- ", nansd(tmp) )
    tmp = float32([clusters[k].dens_string for k=1:k_clust])
    out_df["STRING_DENS"] = nanmean(tmp)
    println( "MEAN STRING DENS: ", out_df["STRING_DENS"][1], " +/- ", nansd(tmp) )
    tmp = float32([clusters[k].meanp_meme for k=1:k_clust])
    out_df["MEME_PVAL"] = nanmean(tmp)
    println( "MEAN MEME LOG10(P-VAL): ", out_df["MEME_PVAL"][1], " +/- ", nansd(tmp) )
    rows = 0; for k in 1:k_clust rows = [rows, clusters[k].rows]; end
    tmp = float32(collect(values(table(rows))))
    out_df["CLUSTS_PER_ROW"] = nanmean(tmp)
    println( "CLUSTS PER ROW: ", out_df["CLUSTS_PER_ROW"][1], " +/- ", nansd(tmp) )
    cols = 0; for k in 1:k_clust cols = [cols, clusters[k].cols]; end
    tmp = float32(collect(values(table(cols))))
    out_df["CLUSTS_PER_COL"] = nanmean(tmp)
    println( "CLUSTS PER COL: ", out_df["CLUSTS_PER_COL"][1], " +/- ", nansd(tmp) )
    out_df
end

function clusters_to_dataFrame( clusters::Dict{Int64,bicluster} )
    out = Array(DataFrame,length(clusters))
    for k in 1:length(clusters)
        b = clusters[k]
        out_r = DataFrame( { "k" => k,
                            "rows" => [join(rownames(ratios)[b.rows],',')],
                            "resid" => b.resid,
                            "dens_string" => b.dens_string,
                            "meanp_meme" => b.meanp_meme,
                            "cols" => [join(colnames(ratios)[b.cols],',')],
                            "meme_out" => [join(b.meme_out,"<<<<>>>>")] 
                            } )
        out[k] = out_r
    end
    rbind(out)
end
