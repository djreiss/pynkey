########## UTILITY FUNCTIONS ################################################################################

## TODO: use HDF5 and JLD to save data (and possibly to use filehash-like feature... see 
##  https://groups.google.com/forum/?fromgroups=#!topic/julia-users/6Kq20HHgOYY and https://github.com/timholy/HDF5.jl
##  (NOTE saving via JLD doesn't work right for my data types yet)

function save_jld(fname, tuple)
    io = GZip.gzopen(fname, "w") ## I get an error right now using gzopen
    try
        serialize(io,tuple) ## doesn't work with gzopen for an array of biclusters so...
    catch
        close(io)
        if endswith(fname, "jldz") ## serialize it uncompressed and then gzip it.
            fname = fname[ 1:(end-1) ]
            io = open(fname, "w")
            serialize(io, tuple)
            close(io)
            run(`gzip -9 -fv $fname`) ## saves about 80% of space  (75% for no -9)
            run(`mv -fv $(fname).gz $(fname)z`)
        end
    end
    close(io)
end

function load_jld(fname)
    io = GZip.gzopen(fname, "r")  ## This can read from the gzipped file as well as the gzopen-ed file
    out = deserialize(io)
    close(io)
    out
end

source = require ## load ## I keep getting hitched up with this so create an alias

const NA = float32( 0 / 0 )

stop() = error( "STOPPING ON PURPOSE" )

## Function to run an external command and get results as a set of strings:
## Oops - just discovered 'readall()' - use that instead
function system( cmd::Cmd )
    ##ll=each_line(cmd) ##`ls`)
    ##lines::Vector{ASCIIString} = []
    ##for i=ll lines=push!(lines,i); end
    lines=readall(cmd)   ## try readsfrom()?
    lines=split(join(lines),'\n')
end

## bsxfun() does this too -- i.e. outer(a,b,+) == bsxfun(+,a,b') -- but actually bsxfun is about 10x slower!
outer{T}(a::Vector{T}, b::Vector{T}, f) = convert( Array{T,2}, [ f(a[i],b[j]) for i=1:length(a), j=1:length(b) ] )

function collect_dict( a, b ) ## Just add them all to a dict, indexed by integer, for now. 
    if ! ( typeof(a) <: Dict ) ## Won't work correctly if each obj. is a dict!
        out = Dict{Int64,typeof(a)}()
        out[1] = a
        a = out
    end
    if ! ( typeof(b) <: Dict ) a[length(a)+1] = b
    else for i=1:length(b) a[length(a)+1] = b[i]; end
    end
    a
end

function table( x )
    d = Dict{typeof(x[1]),Int64}()
    for i=1:length(x) d[x[i]] = get(d,x[i],0) + 1; end
    d
end

##unique( x ) = keys( table( x ) ) ## now in Base

function duplicated( x )
    tab = table( x )
    out = falses( length(x) )
    for i=reverse(1:length(x))
        v = get(tab, x[i], 0)
        if v > 1 
            out[i] = true
            tab[ x[i] ] = v - 1
        end
    end
    out
end

stats{T}(x::Array{T}) = (mean(x),var(x))

function allequal(a,b)
    for i in 1:length(a)
        if a[i] != b[i] return false; end
    end
    return true
end

function remove{T}(a::Array{T},b::T)   ## faster than a[a.!=b] and faster than delete!(a,find(a.==b)), as long as we know b is in a!
    out = Array(T,length(a)-1)
    i = 0
    for x = a
        if x != b out[ i += 1 ] = x; end
    end
    out
end

# function append!{T}(a::Array{T,1}, item::T)
#     ccall(:jl_array_grow_end, Void, (Any, Uint), a, 2)
#     a[end] = item
#     return a
# end

## note [a, [item]] works just as fast as this:
# function append{T}(a::Array{T,1}, item::T)
#     out = Array(T,length(a)+1)
#     i = 0
#     for x = a out[ i += 1 ] = x; end
#     out[end] = item
#     out
# end

import Base.in
function in{T}(a::Array{T}, b::Array{T})  ## in(a,b) same as a%in%b in R
    out = falses(length(a))
    for i = 1:length(a)
        for y = b if a[i] == y out[i] = true; break; end; end
    end
    out
end

# function find_in{T}(a::Array{T}, b::Array{T})  ## faster than find(in(a,b)) ??  NOTE there is already a findin()
#     out::Vector{Int64} = [];
#     broken::Bool = false;
#     for i = 1:length(a)
#         broken = false
#         for y = b
#             if a[i] == y 
#                 out = push(out, i)
#                 broken = true
#                 break
#             end
#             if broken break; end
#         end
#     end
#     out
# end

# function rep{T}(val::T, len::Int)
#     #out = Array(T, len)
#     #fill!(out,val)
#     #out
#     fill(val, len)
# end


## Define which() for vector or matrix that return a 2-column matrix
## Note this is already implemented by find(vector) and findn(matrix)
# function which( x::VecOrMat{Bool} )
#     o=Array(Int,sum(x))
#     j::Int=0
#     for i=1:length(x)
#         if x[i]
#             o[j+=1]=int(i)
#         end
#     end
#     if (length(size(x))==2)
#         o=hcat(map(rem,o,size(x,1)),map(div,o,size(x,1))+1)
#     end
#     o
# end

# ## Index a matrix by a 2-column matrix output by which (above)
# function ref(x::Matrix, inds::Matrix{Int})
#     o=Array(typeof(x[1]),size(inds,1))
#     for i=1:size(inds,1)
#         o[i]=x[inds[i,1],inds[i,2]]
#     end
#     o
# end

## Whew, this is the best way I could find to get rowmeans, rowsds, ignoring NaNs:
## can use amap(f,x,1) to apply f to rows of x but for some reason this isnan() doesn't work with amap():
#isnan{T}(x::Array{T}) = map( isnan, x )

##nanmean(x::Array{Float32}) = mean( x[ ! isnan( x ) ] )

## TODO: should use DataFrames: mean(removeNA(x))  -- but this iterates the vector twice (and 1 copy!)
function nanmean{T}(x::AbstractArray{T}) ## This is a lot faster!
    v=zero(T); nnnan=zero(Int64)
    for xi=x
        if isnan(xi) continue; end
        v += xi
        nnnan += 1;
    end
    v / nnnan
end

function nanvar{T}(x::AbstractArray{T}, m::AbstractArray{T}, corrected::Bool=true)
    n = length(x)
    v=zero(T); xi=zero(T)
    for i=1:n
        xi = x[i]
        if isnan(xi) continue; end
        xi -= m[i]
        v += xi*xi
    end
    return v / (n - (corrected ? one(T) : zero(T))) ##dot(x, x)
end

function nanvar{T}(x::AbstractArray{T}, m::T, corrected::Bool)
    v=zero(T)
    for xi=x
        if isnan(xi) continue; end
        xi -= m
        v += xi*xi
    end
    v / (length(x) - (corrected ? one(T) : zero(T))) ##dot(x, x)
end

nanvar(x::AbstractArray, corrected::Bool=true) = nanvar(x, nanmean(x), corrected)
#nanvar(x::AbstractArray) = nanvar(x, true)

nansd{T}(x::AbstractArray{T}, m::AbstractArray{T}, corrected::Bool=true) = sqrt(nanvar(x, m, corrected))
nansd{T}(x::AbstractArray{T}, m::T, corrected::Bool=true) = sqrt(nanvar(x, m, corrected))
nansd{T}(x::AbstractArray{T}, corrected::Bool=true) = nansd(x, nanmean(x), corrected)
#nansd{T}(x::AbstractArray{T}) = nansd(x, true)
#nansd{T}(x::Vector{T}) = nansd(x, true)

stdize_vector{T}( x::AbstractArray{T} ) =  ( x .- nanmean(x) ) ./ nansd( x ) ## TODO: try Devectorize package???
sdize_vector{T}( x::AbstractArray{T} ) =  x ./ nansd( x )

function range{T}(x::AbstractArray{T}) 
    mn = typemax(eltype(x)); mx = typemin(eltype(x))
    for xi = x
        if xi < mn mn = xi; end
        if xi > mx mx = xi; end
    end
    (mn,mx)
end

##rowmeans(x::AbstractArray) = [ nanmean(x[i,:]) for i=1:size(x,1) ]
function rowmeans{T}(x::AbstractArray{T}) 
    y = Array(T,size(x,1))
    for i=1:size(x,1) y[i] = nanmean(x[i,:]); end
    y
end
##rowsds(x::AbstractArray) = [ nansd(x[i,:]) for i=1:size(x,1) ]
##colmeans(x::AbstractArray) = [ nanmean(x[:,i]) for i=1:size(x,2) ]
##colsds(x::AbstractArray) = [ nansd(x[:,i]) for i=1:size(x,2) ]

function colmeans{T}(x::AbstractArray{T}, filler::Vector{T}) ## pre-allocated vector
    for i=1:length(filler) filler[i] = nanmean(x[:,i]); end
    filler
end

function colmeans{T}(x::AbstractArray{T}) 
    y = Array(T,size(x,2))
    colmeans(x, y)
end

nrow(x::Matrix) = size(x,1)
ncol(x::Matrix) = size(x,2)

## Corrected cor to use equiv of R 'pairwise' (ignoring NAs)
function nancor(x::Vector, y::Vector)
    bad = isnan(x)|isnan(y)
    if (any(bad))
        x=x[!bad]
        y=y[!bad]
    end
    return cor_pearson(x, y)
end

function nancor(x::Matrix)
    o=ones(eltype(x),(size(x,1),size(x,1)))
    for i=1:size(x,1)-1
        xx=vec(x[i,:]) ## apparently this is how you get the output to be a vector
        println(i)
        for j=i+1:size(x,1)
            yy=vec(x[j,:])
            o[i,j]=o[j,i]=nancor(xx,yy)
        end
    end
    return o
end

function reverse_dict( d::Dict ) ## reverse keys <-> elements for a dict
    out = Dict()
    k = keys( d )
    for i=1:length( d ) out[ d[ k[ i ] ] ] = k[ i ]; end
    out
end

function my_sampler( x::Vector, n, probs::Vector )
    ind = 0
    to_samp = [1:length(x)]
    out_rows = zeros(Int, n)
    while ind < n ## add 5 new genes from genes in few clusters -- need to do my own sampler
        tmp_row = Distributions.sample( to_samp )
        if rand() <= probs[tmp_row]
            ind += 1
            out_rows[ind] = x[tmp_row]
        end
        ##println(tmp_row," ",x[tmp_row]," ",ind)
    end
    out_rows
end


## TODO: use pfork from PTools package -- need to write the equiv of a pmap() function.
## "vals" is an input Vector or Dict or otherwise iterable

function pfmap_inner( ind, v, np, f::Function, args... )
    inds = find( mod( [0:(length(v)-1)], np ) .== ind-1 )
    out = [ f(v[ii], args...) for ii=inds ]
    out
end

function pfmap( vals, nprocs, f::Function, args... )
    tmp = pfork( nprocs, pfmap_inner, vals, nprocs, f, args... )
end
