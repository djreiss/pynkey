
type bicluster
    k::Int64
    rows::Vector{Int64}
    cols::Vector{Int64}
    var::Float32
    resid::Float32
    dens_string::Float32
    meanp_meme::Float32
    scores_r::Vector{Float32}
    scores_c::Vector{Float32}
    scores_n::Vector{Float32}
    scores_m::Vector{Float32}
    meme_out::Vector{ASCIIString}
    mast_out::DataFrame ##Array{Any,2}
    changed::BitArray{1} ## rows, cols
end

bicluster( k, rows::Vector, cols::Vector ) = bicluster( k, unique(rows), unique(cols),
                                 ##unique( int64( [rand(1:size(x,2)) for i=1:div(size(x,2),2)] ) ),
                                 typemax(Float32), typemax(Float32), typemax(Float32), typemax(Float32), 
                                 Array(Float32,0), Array(Float32,0), Array(Float32,0), 
                                 Array(Float32,0), Array(ASCIIString,0), DataFrame(), trues(2) )

bicluster( k, rows::Vector, x::NamedMatrix ) = bicluster( k, rows, 
                                 Distributions.sample([1:size(x.x,2)],div(size(x.x,2),2),replace=false,ordered=true) )
