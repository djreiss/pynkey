
## Generate all gapped regexes as pair-wise combo's of kmers, separated by a "gap" of length
## between gap_range[1] and gap_range[2], and possibly with n_ns in each of the k-mers
# function generate_gapped_regexes( kmers::Array{ASCIIString,1}, n_ns, gap_range )
#     if n_ns > 0
#         kmers2 = Dict{ASCIIString,Int64}() ##Array( ASCIIString, length(kmers) )
#         ind = 0
#         for i=1:length(kmers)
#             for k=1:n_ns
#                 for j=1:length( kmers[i] )
#                     km = chars( kmers[i] )
#                     km[j] = 'N'
#                     kmers2[ASCIIString( uint8(km) )] = ind+=1
#                 end
#             end
#         end
#         kmers = [kmers,keys(kmers2)]
#         kmers = [replace(i,"N","[GATC]") for i=kmers]
#     end
#     out = Array(ASCIIString,length(kmers)^2)
#     ind = 0
#     gr1 = gap_range[1]; gr2 = gap_range[2]
#     for i=1:length(kmers)
#         kmi=kmers[i]
#         tmp = "$kmi[GATC]{$gr1,$gr2}"
#         for j=1:length(kmers)
#             kmj=kmers[j]
#             out[ind+=1] = "$tmp$kmj"
#         end
#     end
#     out
# end

## Get hamming distance between 2 strings; short circuit if it becomes more than max.
function hamming_distance( a::ASCIIString, b::ASCIIString, max::Uint8 ) 
    out::Uint8 = 0
    const one::Uint8 = uint8(1)
    for i=1:minimum([length(a),length(b)])
        if a[i] != b[i] out += one; end
        if out > max break; end
    end
    out
end

hamming_distance( a::ASCIIString, b::ASCIIString ) = hamming_distance( a, b, uint8(99) )

require( "extras/sparse.jl" )

## New strategy: for a set of seqs, iterate through all paired k-mers with any gap (between gap_range)
##   add them to a dict. Then group the pairs (i1 i2) (j1 j2) if the hamming distance of (i1,j1)
##   AND of (i2,j2) are BOTH <= n_ns.
function count_gapped_motifs(seqs::Array{ASCIIString,1}, k::Int64, n_ns::Int64,
                             gap_range::Array{Int64,1}, compute_dists::Bool, rc::Bool)
    const max_motif_length = k*2 + gap_range[2]
    d = Dict{ASCIIString,Int64}()
    m1::ASCIIString = ""; m2::ASCIIString = ""; ind::ASCIIString = ""; i2::Int64 = 0; i3::Int64 = 0;
    if rc seqs = [seqs, [ revComp(i) for i=seqs ] ]; end
    for seq in seqs
        for i=1:length(seq)-max_motif_length-1
            if i%5000 == 4999 println(i," ",length(seq)); end
            i2=i+k-1
            m1=seq[i:i2]
            for j=gap_range[1]:gap_range[2]
                i3=i2+j
                m2=seq[i3:i3+k-1]
                ind = "$m1 $m2" ##[m1,m2]
                d[ind] = get(d,ind,0) + 1
            end
        end
    end

    ## group the resulting motifs based on hamming distance described above 
    ## use extras/sparse.jl 
    ## note use of sparse realls slows down the computation (for n_ns>1), 
    ##    so let's use non-sparse and then convert to sparse? OK, but requires RAM!
    ## move hamming into here to short-circuit the computation if distance > n_ns 
    ##    (and then only add it to the distance matrix if it's <= n_ns)
    ## OR: precompute this matrix for ALL similar k-n-k,d motifs and use that matrix?
    ##   do it in test.jl via (motifs,dists)=count_gapped_motifs([seq],4,2,[2,8]);
    if compute_dists
        const n_ns8 = uint8(n_ns) ##; const n_ns81 = uint8(n_ns+1)
        const one::Uint8 = uint8(1); const len::Int64 = k*2+1
        ks = keys(d)
        nk = length(ks)
        println(nk)
        s1::ASCIIString = ""; s2::ASCIIString = ""; dd::Uint8 = 0; ## pre-defining helps?
        dists = zeros(Uint8,nk,nk)  #spzeros
        for i=1:nk-1 
            if i%1000 == 999 println(i," ",nk); end
            s1 = ks[i]
            for j=i+1:nk 
                s2 = ks[j]; dd = 0 ##hamming_distance(s1,s2,n_ns81)
                for k=1:len
                    if s1[k] != s2[k] dd += one; end
                    if dd > n_ns8 break; end
                end
                if dd <= n_ns8 dists[i,j] = dd; end
            end
        end
        dists = sparse(dists)
        gc()
        return( (d,dists) )
    end
    d
end

count_gapped_motifs( seqs ) = count_gapped_motifs( seqs, 4, 1, [2,8], false, true )

## For fun, test the all k-mer generating function
@time for i=1:10 tmp=generate_all_kmers(8); end

## Compute upstream sequences!!!

## Get upstream sequences for genes in cluster 1
seqs=get_sequences(rownames(x)[rows])
kmers4=getBgCounts(seqs[:,3],3)

## First compute ALL k=4,{gap=2-8},k=4(n=2) motifs over genome and their hamming distances 
##    (pairs of 4-mers separated by gap of length 2-8; hamming distances allow up to n=2 
##    mismatches. Do not compute distances as this requires a 65k x 65k matrix!
@time (motifs_all,dists_all)=count_gapped_motifs([seq],4,2,[2,8],false,true);

## OR - do it individually for every sequence!
motifs_all_seqs=Dict{ASCIIString,Any}()
@time begin
for i=1:size(all_seqs,1) 
    println(all_seqs[i,1])
    motifs_all_seqs[all_seqs[i,1]] = count_gapped_motifs([all_seqs[i,3]],4,2,[2,8],true,true); ## Do compute distances
end
end

## Then compute 4{2-8}4(2) motifs for upstream sequences in cluster 1
##   this is pretty fast (including distances) so can compute it on-the-fly
##   otherwise can use dists computed above from entire genome and set param true to false.
## But still need to compute the global counts above to get background.
## Can use this to compute if TOTAL # of motifs in this sequence set is significant
@time (motifs,dists)=count_gapped_motifs(seqs[:,3],4,2,[2,8],true,true);

## Or else use the values computed from each sequence above
## Can use this to compute if # of seqs that have this motif is significant
## e.g.:
best_motif = keys(motifs)[findn(values(motifs).==maximum(values(motifs)))][1]
seq_hits = [ get(motifs_all_seqs[i][1],best_motif,0) for i=seqs[:,1] ]

addprocs_local(3) ## flakey for interactive sessions!
require("junkey/junkey_funcs.jl")
seq_hits_rand = @parallel ( vcat ) for i=1:1000
    println(i)
    samp = randival( 1, size(all_seqs,1), size(seqs,1) )
    ## TO TRY: Instead of recomputing (this next line...), collate values from motifs_all_seqs?
    motifs_rnd=count_gapped_motifs(all_seqs[samp,3],4,2,[2,8],false);
    best_motif_rnd = keys(motifs_rnd)[findmax(values(motifs_rnd))[2]]
    seq_hits_rnd = [ get(motifs_all_seqs[i],best_motif_rnd,0) for i=all_seqs[samp,1] ]
    sum(seq_hits_rnd .> 0)
end
p_value = mean(seq_hits_rand .>= sum(seq_hits .> 0))

