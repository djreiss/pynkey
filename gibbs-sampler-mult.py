##require( "extras/distributions.jl" ) ## for sample()

### ### To run gibbs sampler on a cluster:

### q <- gibbs_site_sampler( seqs=get_sequences( k ), w=20, w_test_range=0,
##        bg_model=mkBgFile( genome_info$all_upstream_seqs, order=0 ), plot="BEST", iter_1=2000, iter_2=20 )

### ### Can re-start (to further optimize) on same motif with arg pssm=q$pssm, or use meme best-fit PSSM via
### ###   arg pssm=meme_scores[[k]]$meme_out[[1]]$pssm
### ### Can find a second motif with arg masks=q$pssm (takes PSSM and masks out seqs best match to that PSSM)
### ### Then display the motif via viewPssm( q$pssm )

### ### "priors" should be a list of 2-by-length(seqs) positional-prior distributions for position of motifs in seqs
### ###   with 2 cols for fwd and rev strand (can be identical of course), and not normalized (e_g_ to weight
### ###   different sequences by different amounts)_
### ### Idea: can increase chance of motifs ONLY in same strand as gene by making col_ 2 of prior lower than col_ 1
### ### Idea: weight different sequences higher if they're more correlated with the bicluster (will only work with
### ###   "zoops" model that will add the max score to the motif counts (rather than 1)) - TODO: is this true?
### ### Idea: use chip-chip as positional prior
### ### Idea: use density of positions of motif in OTHER sequences as prior for given sequence (will require some
### ###   coding intervention within the gibbs_site_sampler function)
### ### All of the above can be done if I change priors from a list to a function (can be dynamically computed
### ###   on the fly)_
### ### TODO: prior on columns in motif as well (and perhaps prior for palindrome?)

## Then for discriminiative motif optimizing:
##  qqq <- optim_pssm( k, "all", pssms=list( q$pssm ) )

#@profile begin

function pssm_to_string( pssm::Array{Float64,2}, cutoff_1, cutoff_2 )
    out = Array( Uint8, size( pssm, 1 ) )
    for i=1:size(pssm,1)
        mx = findmax(pssm[i,:])
        if mx[2] == 0 out[i] = 'n'; continue; end;
        letter = DNA_letters[ mx[2] ]
        value = pssm[ i, mx[2] ]
        if value < cutoff_2 letter = 'n'
        elseif value < cutoff_1 letter = lowercase(letter); end
        out[i] = letter
    end
    ASCIIString( out );    
end

pssm_to_string( pssm ) = pssm_to_string( pssm, 0.7, 0.4 )

function seq_to_res( seq::ASCIIString, rev_comp::ASCIIString )
    res = zeros( Int8, length(seq) )
    for i=1:length(seq) 
        if has( DNA_letter_lookup, seq[i] ) 
            res[i] = DNA_letter_lookup[seq[i]]
        elseif has( iupac_dict, seq[i] ) 
            tmp = iupac_dict[seq[i]]
            res[i] = randi(length(tmp))
            ##if(res[i]>4)println(seq[i]," ",tmp," ",res[i]," ",length(tmp));stop();end
        end
    end
    ##if rev_comp == "no" return( res )
    if rev_comp == "only" res = 5 - reverse( res )
    elseif rev_comp == "both" res = hcat( res, 5 - reverse( res ) ); end
    res
end

seq_to_res( seq ) = seq_to_res( seq, ["no","only","both"][3] )

function subseqs_to_counts( subseqs::Array{Int8,2}, max_scores::Array{Float64,1} ) ##, pseudo ) {
    counts::Array{Float64,2} = zeros( size(subseqs,1), 4 )
    for i=1:size(subseqs,2)
        m = max_scores[ i ];
        for j=1:size(subseqs,1)
            sji = subseqs[j,i]
            if sji > 0 counts[ j, sji ] += m; end ##= counts[ ss( j, i ) - 1 ][ j ] + m;
        end
    end
    counts
end

##function pssm_scan( pssms, scores=NULL, seqs, sums_only=T, pow=0, the_min=1e-30 ) 
##function pssm_scan( pssms, scores, seqs::Array{ASCIIString,1}, res::Array{Array{Int64,2},1}, the_pow::Float64, the_min::Float64, sums_only::Bool ) 
## NOT USED ANYMORE:
# function pssm_scan( pssms::Array{Array{Float64,2}}, scores::Dict{Int64,Array{Float64,2}}, 
#                    seqs::Array{ASCIIString,1}, res::Array{Array{Int64}}, the_pow::Float64, the_min::Float64, 
#                    sums_only::Bool ) 
# #    out_sums::Array{Float64,1} = zeros( length(seqs) ) ## Store the sum of best scores over all PSSMs here
# #    out_maxes::Array{Float64,1} = zeros( length(seqs) ) - 9e9 ## Store the maximum of best scores over all PSSMs here
# #    out_max_posns::Array{Int64,1} = zeros( Int64, length(seqs) ) ## Store position of the the maximum (neg_ for rev_ strand)
#     for pi = 1:length(pssms)
#         const pssm::Array{Float64,2} = pssms[pi]
#         ##const pssmLen::Int64 = size(pssm,1)
#         for si=1:length(seqs)
#             (scores[si][:],max,max_posn) = pssm_scan( pssm, seqs[si], res[si], the_pow, the_min, sums_only ) 
# #            if sums_only == false ## Store the fwd and rev-strand PSSM matches for output
# #                if max > out_maxes[ si ]
# #                    out_maxes[ si ] = max; ## Store the max_ of best scores over all PSSMs
# #                    out_max_posns[ si ] = max_posn; ## and its position
# #                end
# #            end
# #            out_sums[ si ] += max; ## Store the sums of best scores for each seq over all PSSMs
#         end
#     end
# end

function pssm_scan( pssm::Array{Float64,2}, seq::ASCIIString, res::Array{Int64,2}, 
                   in_scores::Array{Float64,2}, sums_only::Bool ) 
    const pssmLen::Int64 = size(pssm,1)
    const sLen::Int64 = length( seq ); ## Its length
    #max::Float64 = -9e9; ## Maximum score over whole seq (over both fwd or rev strands)
    #max_posn::Int64 = 0; ## Position of the maximum (neg_ for reverse strand)
    ##scores::Array{Float64,2} = Array(Float64,sLen,2) ## Store the scores in a 2-column scores matrix
    ssum::Float64 = 1.; sumRev::Float64 = 1.; ind::Int64 = 0; jRange = 1:pssmLen; j::Int64 = 0;
    for i = 0:sLen - pssmLen
        ssum = 1.; sumRev = 1.; ##the_min; sumRev = the_min;
        ##jj = sLen - i + 1
        for j = jRange ##1:pssmLen ##; j ++, jj -- ) {
            ##res::Char = seq[ i + j ]
            ind = res[ i + j, 1 ] ##res == 'G' ? 1 : res == 'A' ? 2 : res == 'T' ? 3 : res == 'C' ? 4 : 0
            ##if ind > 0 
            ssum *= pssm[ j, ind ]##; end; ## Compute PSSM match score at posn i (fwd strand)
            ##res = seq[ jj -= 1 ]
            ##ind = res[ jj -= 1, 2 ] ##res == 'G' ? 4 : res == 'A' ? 3 : res == 'T' ? 2 : res == 'C' ? 1 : 0
            ind = res[ i + j, 2 ] ##res == 'G' ? 4 : res == 'A' ? 3 : res == 'T' ? 2 : res == 'C' ? 1 : 0
            ##if ind > 0 
            sumRev *= pssm[ j, ind ]##; end; ## Compute PSSM match score at posn i (rev strand)
        end
        i += 1;
        if sums_only == false ## Store the PSSM match scores
            in_scores[ i, 1 ] = ssum
            in_scores[ i, 2 ] = sumRev
        end
        #if sum > max
        #    max = sum; max_posn = i; ## Store the maximum seen over this sequence
        #end
        #if sumRev > max
        #    max = sumRev; max_posn = i; ##-i - 1; ## ___either fwd or reverse
        #end
    end
    #(scores, max, max_posn)
    in_scores
end

function gibbs_site_sampler( seqs, w, model, iter_1, iter_2, n_start, n_retries, do_sample, pseudo,
                            bg_model, bg_order, bg_pow, w_range, w_select, quiet, pssm, masks, priors, parallel )
    ##seqs_rc <- rev_comp( seqs )
    if typeof(pssm) != Nothing ## user input a starting "seed" pssm or motif string
        if typeof( pssm ) == ASCIIString
            str = pssm
            pssm = zeroes( length( str ), 4 )
            for i=1:length(str) pssm[i, DNA_letter_lookup[ uppercase( stri[i] ) ]] = 1; end
        end
        if any( pssm .== 0 )
            pssm += pseudo
            for i=1:size(pssm,1) pssm[i,:] /= sum(pssm[i,:]); end
        end

        if w != size( pssm, 1 )
            ##if w < size( pssm, 1 ) pssm = pssm[ 1:w, : ]
            ##elseif w > size( pssm, 1 ) pssm = hcat( pssm, zeroes( w-size( pssm, 1 ), 4 ) + 0.25 ); end
            w = size( pssm, 1 )
        end
    end

    n_seqs = length( seqs )
    seq_lens = int32([ length( seqs[ i ] ) for i=1:n_seqs ])
    max_len = max( seq_lens )
    w_range += w ## w_range is the offset from the input w that can be altered.

    ##if ( is_null( bg_model ) ) 
    bg_freqs = nothing
    if bg_order >= 0
        bgCounts = getBgCounts( seqs, [0:bg_order], false ); ## 0th order
        ##bgCounts = del( del( bgCounts, "S" ), "Y" ); ## TODO: handle N/S/Y/other codes more gracefully
        bgFreqs = getBgFreqs( bgCounts );
        bg_model = bgFreqs; bg_model[ "order" ] = bg_order;

        ## Get freqs along each sequence for bg model
        ##bg_freqs = NULL ##list()
        ##if ( ! is_null( bg_model ) ) {
        bg_freqs = Dict{Int64,Array{Float64,2}}() ##list() ##bg_freqs_smooth = list()
        ##bg_model = unlist( bg_model, use_names=T )
        for ord = 0:bg_order
            for s = 1:n_seqs
                seq = seqs[ s ]
                rseq = revComp(seq)
                nc = seq_lens[ s ]
                tmp1 = float64( [ get(bg_model, seq[i:i+ord], 1) for i=1:nc-ord ] )
                tmp2 = float64( [ get(bg_model, rseq[i:i+ord], 1) for i=1:nc-ord ] )
                tmp = hcat( tmp1, tmp2 )
                ##tmp = cbind( "for"=bg_model[ substring( seq, 1:nc, ( order + 1 ):( nc + order ) ) ],
                ##           bg_model[ substring( rev_comp( seq ), 1:nc, ( order + 1 ):( nc + order ) ) ] )
                ##tmp[ is_na( tmp ) ] = 1
                if ~has( bg_freqs, s ) bg_freqs[ s ] = tmp;
                else bg_freqs[ s ][ 1:size(tmp,1), : ] .*= tmp; end
                ##bg_freqs[s][:] = 1.0
            end
        end

        for s = 1:n_seqs
            tmp = bg_freqs[ s ]
            for i=1:2 tmp[:,i] = (tmp[:,i] ./ max(tmp[:,i])).^bg_pow; end
        end
    end

    best_MAP = -9e9
    start_pssm = pssm
    best_out = pssm
    start_w = w

#    scores = Dict{Int64,Array{Float64,2}}()
#    for i = 1:n_seqs scores[ i ] = zeros( seq_lens[i], 2 ); end

    res = [ seq_to_res( seq, "both" ) for seq=seqs ];
    ##return({"seq"=>seqs[1],"res"=>res[1],"bg_freqs"=>bg_freqs[1],"bgFreqs"=>bgFreqs,"bgCounts"=>bgCounts})

  ## Mask out best locations of "mask" motifs in each seq
#   if ( ! is_null( masks ) ) {
#     if ( is_matrix( masks ) ) masks = list( masks )
#     for ( mask in masks ) {
#       ## Get the log2-match scores at each location in each sequence to the current pssm
#       if ( any( mask == 0 ) ) {
#         mask = mask + pseudo ##/ length( seqs )
#         mask = scale( mask, colSums( mask ), center=F )
#         ##mask = apply( mask, 2, function(i) i / sum( i ) )
#       }
#       ##scores = ##pssm_match_scores( res, mask, rev="both", pseudo=0, logify=T )
#       pssm_scan_cpp( list( log2( mask ) ), scores, seqs, sums_only=F, pow=2, min=1e-30 ) ##$scores

#       ww = nrow( mask )
#       for ( i in 1:n_seqs ) {
#         max_len = seq_lens[ i ] - ww + 1
#         scores_2 = scores[[ i ]][ 1:max_len, ] ##cbind( "for"=scores[[ i * 2 - 1 ]], scores[[ i * 2 ]] )
#         ## Convert to a freq, and divide by the background freqs
#         if ( ! is_null( bg_freqs ) ) scores_2 = scores[[ i ]] / bg_freqs[[ i ]][ 1:max_len, ]
#         ## Get the max freq between for- and rev- strand scores for each site
#         posn = which_max( scores_2 ) ## Get position of highest scores in each seq
#         strand = if ( posn > nrow( scores_2 ) ) 2 else 1 ##which_max( scores[[ i ]][ ,posn ] )
#         ##if ( strand == 2 ) posn = posn - nrow( scores_2 )
#         if ( strand == 2 ) posn = nrow( scores_2 ) - ( posn - nrow( scores_2 ) ) + 1 ##posn - nrow( scores_2 )
#         start = max( 1, posn - w )
#         end = min( posn + ww, max_len ) ##nrow( bg_freqs[[ i ]] ) )
#         ## Do it by setting the bg_freq for these posns to a really big number
#         if ( ! is_null( bg_freqs ) ) bg_freqs[[ i ]][ start:end, ] = Inf
#       }      
#     }
#   }

    ##max_scores = rep( 0, n_seqs )
    ##for ( iter1 in 1:iter_1 ) {
    ##qqq = apply_func( 1:iter_1, function( iter1, ___ ) {
    ##counts = nothing ## for debugging, make it global

    max_scores = zeros( n_seqs )
    posns = zeros( Int16, n_seqs ) ##as_integer( runif( n_seqs, 1, max_lens - 1 ) )
    seq_scores = [ Array(Float64,length(seq),2) for seq=seqs ]
    
    for iter1=1:iter_1
        ##if ! quiet println("ITER1: ",iter1); end
        if ! quiet print(iter1,"..."); end #flush(stdout_stream); end
        best_MAP_2 = -9e9
        w = start_w
        max_lens = seq_lens - w + 1
        subseqs = zeros( Int8, w, n_seqs )
        inds = 0:( w - 1 )
        tries = 1
        shft = wind = 0
        max_scores[:] = 0
        posns[:] = 0

      ## Start with random PSSM from randomly-sampled sites on n_start random seqs (or sampled from "priors")
        if typeof(start_pssm) != Nothing ##|| w != w ## w==w ==> !is.na(w)
            pssm = similar(start_pssm)
            for i=1:size(pssm,1) pssm[i,:] = start_pssm[i,:] / sum(start_pssm[i,:]); end
        else
        ## Make sure chosen start sites are not "masked"
            for i=randi((1,n_seqs),n_start)
                ##bgf = if ( is_null( bg_freqs ) || is_null( bg_freqs[[ i ]] ) ) t( 1 )
                ##    else bg_freqs[[ i ]][ 1:( nchar( seqs[ i ] ) - w - 1 ), ]
                prob = 1
                if priors != nothing 
                    prob = priors( i, posns, seqs, w )[ :,1 ] ##; prob /= sum(prob)
                end
                if bg_freqs != nothing && has(bg_freqs,i) 
                    bgf = bg_freqs[ i ][ 1:( length( seqs[ i ] ) - w - 1 ), : ]
                    prob ./= bgf
                end
                posn = sample( [1:length( prob )], prob[:] ) ##posns[ i ]
                strand = posn > size( bgf, 1 ) ? 2 : 1
                ##if strand == 2 posn = size( bgf, 1 ) - ( posn - size( bgf, 1 ) ) + 1; end ##posn - nrow( scores_2 )
                if strand == 2 posn -= size( bgf, 1 ); end
                ## Get subseqs same length as motif width, at each site in each seq
                subseqs[ :,i ] = res[ i ][ posn + inds, strand ]
                max_scores[ i ] = 1
                ##println(subseqs)
                ##sampled_posns = rbind( sampled_posns, c( iter1, i, posn, strand ) )
            end  
        
            subseqs[ isnan(subseqs) ] = randi( (1,4), sum( isnan(subseqs) ) )
            ## Count up # of occurences of each residue in each position in the subseqs; add pseudocounts
            counts = subseqs_to_counts( subseqs, max_scores ) + pseudo
            ## Convert to a frequency
            pssm = similar(counts)
            for i=1:size(pssm,1) pssm[i,:] = counts[i,:] / sum(counts[i,:]); end
        end ## if

        iter2 = 0
        posns = zeros( Int16, n_seqs ) ##as_integer( runif( n_seqs, 1, max_lens - 1 ) )

        ##for iter2=1:iter_2
        while iter2 <= iter_2
            iter2 += 1
            if wind != 0
                w = w + wind
                if w < w_range[ 1 ] w = w_range[ 1 ]
                elseif w > w_range[ 2 ] w = w_range[ 2 ]; end
                max_lens = seq_lens - w + 1
                subseqs = zeros( Int8, w, n_seqs )
                inds = 0:( w - 1 )
                posns[:] = 0
            end

            ml = seq_lens - size( pssm, 1 ) + 1

            for i=1:n_seqs 
                scores = pssm_scan( pssm, seqs[i], res[i], seq_scores[i], false )
                scores_2 = scores[ 1:ml[ i ], : ]
                ## Adjust scores by the background freqs
                if bg_freqs != nothing
                    bgf = bg_freqs[ i ]
                    if size( scores_2, 1 ) > size( bgf, 1 ) scores_2 = scores_2[ 1:size( bgf, 1 ), : ]
                    elseif size( scores_2, 1 ) < size( bgf, 1 ) bgf = bgf[ 1:size( scores_2, 1 ), : ]; end
                    scores_2 ./= bgf
                end
                ## Adjust scores by priors, if available
                if priors != nothing
                    prior = priors( i, posns, seqs, w ) 
                    if size( scores_2, 1 ) < size( prior, 1 ) prior = prior[ 1:size( scores_2, 1 ), : ]
                    elseif size( scores_2, 1 ) > size( prior, 1 ) scores_2 = scores_2[ 1:size( prior, 1 ), : ]; end
                    scores_2 .*= prior
                end

                ## Get position of highest scores in each seq by sampling (or max if desired)
                if do_sample
                    ##scores_2 /= sum(scores_2)
                    posn = sample( [1:length( scores_2 )], scores_2[:] )
                else ## faster if not sampling
                    posn = findmax( scores_2 )[ 2 ]
                end
                ## Get the strand that each position's score corresponds to
                strand = posn > size( scores_2, 1 ) ? 2 : 1
                ## If zoops, the counts will be scaled by the proportion of the prob. that is in this position
                if model == "zoops" ## Weight it by the sum of the prior (allow different weights for diff. seqs)
                    ms = scores_2[ posn ] / sum( scores_2 )
                    if ( priors != nothing ) ##&& ! is.null( priors[[ i ]] ) ) ms <- ms * sum( priors[[ i ]] )
                        ms .*= sum( priors( i, posns, seqs, w ) ); end
                    max_scores[ i ] = ms
                else
                    max_scores[ i ] = 1
                end
                ## Save the positions (neg_ if strand is rev_ comp_)
                if strand == 2 posn -= size(scores_2, 1); end ##posn = size( scores_2, 1 ) - ( posn - size( scores_2, 1 ) ) + 1; end
                posns[ i ] = strand == 1 ? posn : -posn
                ## Get subseqs same length as motif width, at each site in each seq, shifted if shift != 0
                if shft != 0 posn += shft; end
                t_inds = posn + inds
                if shft != 0 || wind != 0 ## Check for out-of-bounds problems with shifting or resizing
                    if t_inds[ 1 ] < 1 
                        t_inds = t_inds - t_inds[ 1 ] + 1
                    elseif t_inds[ end ] > size( res[ i ], 1 ) 
                        t_inds = t_inds - ( t_inds[ end ] - size( res[ i ], 1 ) )
                    end
                end
                subseqs[ :,i ] = res[ i ][ t_inds, strand ]
            end ## iter i over seqs

            ## Count up # of occurences of each residue in each position in the subseqs; add pseudocounts
            counts = subseqs_to_counts( subseqs, max_scores ) + pseudo

            ## Convert to a frequency
            pssm = similar(counts)
            for i=1:size(pssm,1) pssm[i,:] = counts[i,:] / sum(counts[i,:]); end

            ## Compute new MAP score to be maximized
            if w_select == "MAP" 
                MAP = mean( pssm .* counts )
            elseif w_select == "BIC"
                BIC = -2 * sum( pssm .* counts ) + 3 * size( pssm, 1 ) * log( sum( counts[ 1,: ] ) )
                MAP = -BIC
            elseif w_select == "AIC"
                AIC = -2 * sum( pssm .* counts ) + 2 * 3 * size( pssm, 1 )
                MAP = -AIC
            end

            ##pssm = new_pssm ##copy( new_pssm )

            if iter1 == 1 && iter2 == 1 
                best_out = {"pssm"=>pssm, "counts"=>counts-pseudo, "iters"=>[0,iter2,tries],
                            "posns"=>hcat( posns, max_scores ), "MAP"=>MAP}
            end

            ## Save pssm if MAP is best ever seen
            if MAP > best_MAP #* 1.01 ##00001
                if ! quiet 
                    if best_out["iters"][1] != iter1 println(); end
                    @printf( "*** %4d %4d %4d %4d %4d %4d\tMAP = %.3f\t%.3f\t%.3f\t%s\n", iter1, iter2,
                            tries, wind, shft, size( pssm, 1 ), best_MAP, MAP, best_MAP_2, pssm_to_string( pssm ) )
                end
                best_MAP = MAP
                best_out = {"pssm"=>pssm, "counts"=>counts-pseudo, "iters"=>[iter1,iter2,tries],
                         "posns"=>hcat( posns, max_scores ), "MAP"=>MAP}
            end
            if MAP > best_MAP_2 ## got a better motif (local iter) so keep trying
                best_MAP_2 = MAP
                tries = 1
                shift = 0
                iter2 -= 1
            elseif tries < n_retries ## didn't get a better motif (local iter) so try shift/resize until we get a better one
                tries += 1
                samp = [- 1, 0, 1 ]
                if tries > n_retries / 2
                    if shft == 0 shft = sample( samp )
                    else shft = sample( [ shft, shft, samp ] ); end
                    ##samp[ samp != shft ], 1 )
                    if diff( w_range )[1] > 0 && wind == 0
                        wind = sample( samp ); ttries = 0
                        while ttries < 10 && ( w + wind < w_range[ 1 ] || w + wind >= w_range[ 2 ] )
                            wind = sample( samp )
                            ttries = ttries + 1
                        end
                        if ttries >= 10 wind = 0; end
                    ##end ##else wind = 0
                    else 
                        wind = sample( [ wind, wind, samp ] ); 
                    end
                    ##samp[ samp != wind ], 1 )
                else 
                    wind = shft = 0
                end
            else ## didn't get a better motif (local iter) and already tried shifting/resizing n_retries times, so break out.
                break
            end
#            if MAP >= best_MAP_2 #* 1.01 ##00001
#                best_MAP_2 = MAP
#                tries = 1
#                shift = 0
#            end
        end
    end
    sites = Array( ASCIIString, length(seqs) )
    w_final = size(best_out["pssm"],1)
    for s in 1:length(seqs)
        posn = best_out["posns"][s,1]
        sites[s] = posn > 0 ? seqs[s][posn:posn+w_final-1] : revComp(seqs[s])[abs(posn):abs(posn)+w_final-1]
    end
    best_out["sites"] = sites''
    ##best_out["posns"] = hcat( best_out["posns"], sites'' )
    best_out
end

## Note, higher pseudocount works better with model="zoops" (i.e. to ident. which seqs don't have the motif)
gibbs_site_sampler( seqs, pssm ) = gibbs_site_sampler( seqs, 18, ["zoops","oops"][1], 100, 50, 1, 20, true,
                                                      0.1, ## 1/length(seqs),
                                                      nothing, 3, 1.0, [-4,4],
                                                      ["BIC","MAP","AIC"][ 1 ], false, 
                                                      pssm, nothing, nothing, true )

gibbs_site_sampler( seqs ) = gibbs_site_sampler( seqs, nothing )

#end # profile

# w=18; model="zoops"; iter_1=100; iter_2=50; n_start=1; n_retries=20; do_sample=true;
# pseudo=0.01 ##length(seqs)/100; 
# bg_model=nothing; bg_order=0; bg_pow=1.0; w_range=[-4,4]+w;
# w_select=["BIC","MAP","AIC"][ 1 ]; quiet=false; ##plot=c("BEST";"TRUE";"LAST";"FALSE");
# pssm=nothing; masks=nothing; priors=nothing; parallel=true ##)

gibbs_test = function(genes)##,anno,genome_seq,op_table)
    seqs = get_sequences(genes)[:,3]
    seqs = keys( table( seqs ) ) ## uniquify the seqs

    #srand(10)

    if false
        insertion = "GAATACGATGAG"
        for i=1:length(seqs)
            seq = seqs[i]
            pos = 50 ##sample([1:length(seq)-12])
            rc = sample([-1,0,0,1])
            println(pos," ",rc)
            if rc == 0 continue; end
            insrt = rc == 1 ? insertion : revComp(insertion)
            s1 = seq[1:pos-1]; s2 = seq[pos+11:end]
            seqs[i] = "$s1$insrt$s2"
        end
    end
end

#@time gibbs_out = gibbs_site_sampler(seqs)
#@time gibbs_out2 = gibbs_site_sampler(seqs, gibbs_out["pssm"])
