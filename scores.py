## Trial: 1.0*r - 0.9*n + 5.0*m + 0.3*v ... NOTE n has NOT been corrected so we need to negate it here!
## DONE: add volume v0 factor (up-weight to grow clusters)
## DONE: add gene popularity g0 factor (down-weight for adding genes that are already in lots of other clusters and
##      up-weight for genes that are in very few)        
## DONE: do we need weight_c? Might be nice to upweight gains for columns vs. rows or vice versa.
## TODO: Load and use other networks (with their respective scores)

## Note can make these weights a global var and update them in this function.
## Just need to define them as global IN the function
function get_score_weights(iter)
    global n_iters, ratios, max_network_weight, max_motif_weight
    const mn = max_network_weight
    const mm = max_motif_weight
    weight_r =  1.0
    weight_n =  0.0 + -mn * float32(iter-1) / n_iters   ## increase linearly from 0 at iter=1 to 0.9
    weight_m =  1.0 +  mm * float32(iter-1) / n_iters * ((iter<=5) ? 0 : 1) ## ramp up from 1 to 1.8 starting at iter=6
    weight_c =  0.0 + 0.2 * size(ratios.x,1)/size(ratios.x,2)/12.0 ## ??? ## 1.2 works good for hpy
    weight_v =  0.1 + 0.3 * float32(iter-1) / n_iters  ## ramp up from 0.3 to 0.8
    weight_g =  0.1 + 0.1 * float32(iter-1) / n_iters  ## ramp up from 0.3 to 0.8
    (weight_r, weight_n, weight_m, weight_c, weight_v, weight_g)
end

function get_score_weights()
    global iter
    get_score_weights(iter)
end

function get_combined_scores( r::Vector{Float32}, n::Vector{Float32}, m::Vector{Float32}, v::Vector{Float32}, g::Vector{Float32} )
    ##(weight_r, weight_n, weight_m, weight_c, weight_v, weight_g) = get_score_weights()
    out::Vector{Float32} = Array(Float32, length(r))
    for i in 1:length(out) 
        out[i] = get_combined_score( r[i], n[i], m[i], v[i], g[i] )
#                 weight_r * ( isnan(r[i]) || isinf(r[i]) ? 0 : r[i] ) +  ## twice as fast to unroll it here
#                 weight_n * ( isnan(n[i]) || isinf(n[i]) ? 0 : n[i] ) +  ##  but it is already really fast so
#                 weight_m * ( isnan(m[i]) || isinf(m[i]) ? 0 : m[i] ) +  ## do it the simpler way
#                 weight_v * ( isnan(v[i]) || isinf(v[i]) ? 0 : v[i] ) +
#                 weight_g * ( isnan(g[i]) || isinf(g[i]) ? 0 : g[i] )
    end
    out
end

function get_combined_score( r::Float32, n::Float32, m::Float32, v::Float32, g::Float32 )
    (weight_r, weight_n, weight_m, weight_c, weight_v, weight_g) = get_score_weights()
    weight_r * ( isnan(r) || isinf(r) ? 0 : r ) + 
    weight_n * ( isnan(n) || isinf(n) ? 0 : n ) + 
    weight_m * ( isnan(m) || isinf(m) ? 0 : m ) +
    weight_v * ( isnan(v) || isinf(v) ? 0 : v ) +
    weight_g * ( isnan(g) || isinf(g) ? 0 : g )
end

## Up-weight moves into clusters with <= 15 rows; down-weight moves out of clusters with >15 rows
## DONE: a little less heuristic?
## DONE? less "sharply-peaked" at 15? Perhaps make it "level out" between say 8 and 23
## DONE: plot this and see -- I think it DECREASES for really big volumes, want to fix this
function get_cluster_volume_row_scores( b::bicluster, is_in_r )
    global avg_genes_per_cluster
    const thresh = avg_genes_per_cluster
    lr = length(b.rows)
    score_vr = float32( [ ( is_in_r[i] ? +1 : -1) * ( thresh - lr ) for i in 1:length(is_in_r) ] )
    if lr >= thresh - 7 & lr <= thresh + 7 score_vr /= 5
    elseif lr >= thresh - 12 & lr <= thresh + 12 score_vr /= 2
    end
    score_vr
end

## Just up-weight moves that add columns (to prevent shrinkage)
## NOTE that as is, the weighting decreases for really high #cols... is this what we want?
function get_cluster_volume_col_scores( b::bicluster, is_in_c )
    lc = length(b.cols)
    score_vc = float32( [ (is_in_c[i] ? -1 : +1) / lc for i in 1:length(is_in_c) ] )
    ##score_vc = fill(NA, length(b.scores_c)) ## Do we really care how many cols there are?
    score_vc
end

## Scores for number of clusters each gene is in
function get_cluster_row_count_scores( counts_g::Vector{Int32} )
    global avg_clusters_per_gene
    const thresh = avg_clusters_per_gene ##1.3 ## 2.0 ## 3.0 ## lower is better; coerce removing if gene is in more than 2 clusters
    float32(counts_g - thresh)
end

## Up-weight moves OUT if counts_g is HIGH, and moves IN if counts_g is LOW
function get_cluster_row_count_scores( cluster::bicluster, counts_g::Vector{Int32}, is_in_r::BitArray )
    score_g = get_cluster_row_count_scores( counts_g )
    score_g = float32( [ (is_in_r[i] ? -1 : +1) * score_g[i] for i in 1:length(is_in_r) ] ) 
    score_g
end

function get_cluster_row_count_scores( cluster::bicluster, counts_g::Vector{Int32} )
    inds = [1:length(cluster.scores_r)]
    is_in_r = in( inds, cluster.rows )
    get_cluster_row_count_scores( cluster, counts_g, is_in_r )
end

##get_cluster_expr_rowcol_scores{T}( b::bicluster, x::NamedMatrix{T} ) = get_cluster_expr_rowcol_scores( b, x, false )

function get_cluster_expr_rowcol_scores{T}( b::bicluster, x::NamedMatrix{T} ) ##, verbose )
    const var_add = 0.1
    ## Try getting the effect of adding/removing each row/col from this cluster on its total variance
    ## DONE: want to have biclusters with larger row-variance, so normalize var(x)/var(mn)
    ## Other variables to be tweaked: currently using out_var(x)=var(x)/(var(colmeans(x))+1) --
    ##       perhaps the "1" should be lower... see "var_add" variable
    ## Commented out volume factor - shoul add this in to total, but externally to this function
    rows=b.rows; cols=b.cols; score_r=b.scores_r; score_c=b.scores_c

    ##const bsize = length(rows)*length(cols)
    const vr = bicluster_var( b, x )
    resid = bicluster_residue( b, x )

    ##xx::Matrix{T} = x.x[:,cols] ## This uses less RAM than previous version; about same speed
    xx::ArrayView{Float32,2,Array{Float32,2},(Array{Int64,1},Array{Int64,1})}=view(x.x,[1:size(x.x,1)],cols) ## This uses less RAM; perhaps a tiny bit slower
    mn::Vector{T} = colmeans(xx[rows,:])
    x2::Matrix{T} = similar(xx) ## Get the full matrix; subtract the mean biclust profile
    for i=1:size(x2,2) x2[:,i] = xx[:,i] - mn[i]; end
    
    ## Get changes in variance for adding/removing rows; normalize by change
    ## in cluster volume. Lower is better!
    ##v_factor::T = length(cols) / bsize
    if length(score_r) == 0 score_r = Array(T,size(x.x,1)); end
    tmpVec::Vector{T} = Array(T, size(xx,2))
    for r in 1:size(x.x,1) ## Iterate over rows
       isIn = in(rows, r) ##any(rows .== r) ## r is in the bicluster
       newR = isIn ? remove(rows,r) : [rows, [r]] ##append(rows, r) ##[rows,r] ## CANT DO: prevent reallocation of vector here?
       newvr = nanvar(x2[newR,:]) / ( nanvar(colmeans(xx[newR,:], tmpVec)) + var_add ) ## DONE -- prevent realloc of matrix here
       ## NOTE: the colmeans() above is the slow part here -- convince yourself that you need it!
       score_r[r] = ( newvr - vr ) #/ (vr+0.01) ##+ ( isIn ? v_factor : -v_factor )
    end

    ##xx=x.x[rows,:]
    xx=view(x.x,rows,[1:size(x.x,2)])
    mn=colmeans(xx)
    x2=similar(xx) ## Get the full matrix; subtract the mean biclust profile
    for i=1:size(x2,2) x2[:,i] = xx[:,i] - mn[i]; end

    ##v_factor = length(rows) / bsize
    if length(score_c) == 0 score_c = Array(T,size(x.x,2)); end
    for c in 1:size(x.x,2) ## Iterate over cols
       isIn = in(cols, c) ##any(cols .== c) ## c is in the bicluster
       newC = isIn ? remove(cols,c) : [cols, [c]] ##append(cols, c) ##[cols,c]
       newvr = nanvar(x2[:,newC]) / ( nanvar(mn[newC]) + var_add )
       score_c[c] = ( newvr - vr ) #/ (vr+0.01) ##+ ( isIn ? v_factor : -v_factor )
    end

    b.var = vr
    b.resid = resid
    b.scores_r = score_r ## Not needed since the vector is filled (passed by reference)
    b.scores_c = score_c ## but we'll do it here anyway b/c it doesn't slow things down much
    b
end

#@iprofile begin
## NOTE WE ASSUME THE NETWORK IS "SYMMETRIZED" if it is undirected - i.e. (n1,n2,w) and (n2,n1,2) are both in the dataframe.
## NOTE we return log10(density) changes, which is OPPOSITE of ROW/MOTIF scores because increases are better
function get_cluster_network_row_scores( b::bicluster, network::DataFrame ) 
    global ratios, all_genes
    ## Note that right now rows in the bicluster is the index into the rownames of the expression matrix!!! TODO: change this!
    r_rownames = rownames(ratios)
    rows = r_rownames[b.rows] ##keys( all_genes )[ in( values( all_genes ), b.rows ) ]
    net = sub( network, findin( network["x1"], rows ) ) ## Assume network is symmetric!
    grps = groupby( net, "x2" ) ## each group is a subnetwork with all genes in bicluster (x1) connected to a given gene
    grpnames = convert( Vector{ASCIIString}, [ grps[i]["x2"][1] for i=1:length( grps ) ] ) ## the given gene for each group
    grpname_lookup = [ grpnames[i] => i for i=1:length(grps) ]
    new_net = net2 = sub( net, findin( net["x2"], rows ) )   ## subnetwork just connecting between nodes in the bicluster
    sum_weights = sum( new_net["x3"] )
    dens = log10( sum_weights / length(rows)^2 + 1e-9 ) ## Already symmetrized, need to decrease count by 1/2
    n2 = net2["x2"] ##Vector{ASCIIString} = net2[2]

    score_n = b.scores_n
    if length(score_n) == 0 score_n = Array(Float32,size(ratios.x,1)); end    
    new_dens = 0.
    for r in r_rownames
        isIn = in(rows, r) 
        if isIn ## r is in the bicluster -- remove the node from the bicluster-only subnetwork and recalc the density
            newR = remove(rows,r)
            new_net = sub( net2, findin( n2, newR ) ) ##net2[2].data, newR ) ) )
            new_dens = sum( new_net["x3"] ) / length(newR) / length(rows)
        else
            try 
                new_net = grps[ grpname_lookup[r] ]
                new_dens = ( sum_weights + sum( new_net["x3"] ) ) / (length(rows)+1) / length(rows)
            catch ## no subnetwork for this "r"; so use the bicluster's subnetwork, but diminished by increased length(newR)
                new_dens = sum_weights / (length(rows)+1) / length(rows)
            end
        end
        score_n[all_genes[r]] = log10(new_dens+1e-9) - dens ##/ (dens+0.1) ##+ ( isIn ? v_factor : -v_factor )  ## lower is BETTER
    end
    b.scores_n = score_n
    b.dens_string = dens
    b
end
#end # profiler

function get_cluster_meme_row_scores( b::bicluster )
    global ratios, all_genes
    #println("IN HERE: MOT SCORES $(b.k)")

    score_m = b.scores_m
    if length(score_m) == 0 score_m = Array(Float32,size(ratios.x,1)); end
    r_rownames = rownames(ratios)
    genes = r_rownames[b.rows]
    #sz=size(b.mast_out);#println("HERE MOT1 $sz")
    ## TODO: don't only use values for genes that are in the ratios matrix.
    ## DONE: make mast_out into a DataFrame for easier searching, subsetting
    if size(b.mast_out,1) <= 0
        b.scores_m = float32( zeros( size(ratios.x,1) ) )
        b.meanp_meme = NA
        return( b )
    end
    df = sub( b.mast_out, findin( b.mast_out["Gene"], genes ) )
    #sz=size(df);#println("HERE MOT2 $sz")
    mn::Float32 = nanmean( log10( df["P-value"].data ) )
    pvals::Dict{ASCIIString,Float32} = Dict{ASCIIString,Float32}()
    for i in 1:size(b.mast_out,1) pvals[b.mast_out["Gene"].data[i]] = b.mast_out["P-value"].data[i]; end
    not_there = r_rownames[ ! in(r_rownames, b.mast_out["Gene"].data) ]
    for r in not_there pvals[r] = NA; end
    for g in r_rownames ## Iterate over rows
       isIn = in(genes, g) ##any(rows .== r) ## r is in the bicluster
       newR = isIn ? remove(genes,g) : [genes,[g]]
       pvs = float32( [ pvals[r] for r in newR ] )
       # pvs = df["P-value"].data
       # if isIn pvs = sub( df, df["Gene"].data .!= g )["P-value"].data
       # else    pvs = [ pvs, sub( b.mast_out, b.mast_out["Gene"].data .== g )["P-value"].data ]
       # end
       newmn = nanmean( log10( pvs ) )
       score_m[all_genes[g]] = newmn - mn ##/ (mn+0.1) ##+ ( isIn ? v_factor : -v_factor )
    end

    b.scores_m = score_m
    b.meanp_meme = mn
    b
end

