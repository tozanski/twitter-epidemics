module TomekUtils

using DataStructures
using SparseArrays

export compress32
export countunique
export countuniquesorted
export arrecdf
export trapzdx!
export trapzdx
export float_bisection
export find_crit

compress32(arr::SparseMatrixCSC{T,I} ) where {T, I<:Integer} = 
    SparseMatrixCSC{UInt8, UInt32}(
        arr.m |> UInt32,
        arr.n |> UInt32,
        arr.colptr .|> UInt32,
        arr.rowval .|> UInt32,
        arr.nzval)

function countunique(arr::AbstractArray{T,1}) where T
    d = Dict{T,Int}()
    for val in arr
        count = get(d, val, 0)
        d[val] = count + 1
    end
    
    N = length(d)
    vals = Array{T,1}()
    sizehint!(vals,N)
    counts = Array{Int,1}()
    sizehint!(counts,N)
    for (val,count) in d
        push!(vals,val)
        push!(counts,count)
    end
    
    return vals,counts
end


function countuniquesorted(arr::AbstractArray{T,1}) where T
    d = Dict{T,Int}()
    for val in arr
        count = get(d, val, 0)
        d[val] = count + 1
    end
    
    N = length(d)
    vals = Array{T,1}()
    sizehint!(vals,N)
    counts = Array{Int,1}()
    sizehint!(counts,N)
    for (val,count) in sort(collect(d))
        push!(vals,val)
        push!(counts,count)
    end
    
    return vals,counts
end

function arrecdf(data::AbstractArray{V,1}) where V<: Number
  x, counts = countuniquesorted(data)
  cumcounts = cumsum(counts)
  p = cumcounts // cumcounts[end]
  return x, p
end

function arrepdf(data::AbstractArray{V,1}) where V <: Number
  x, counts = countuniquesorted(data)
  p = counts // sum(counts)
  return x, p
end


function empdist2D(W1::AbstractVector{V1}, W2::AbstractVector{V2}) where {V1<:Real,V2<:Real}
    W01 = unique(W1)
    W02 = unique(W2)
    sort!(W01)
    sort!(W02)
    
    I = Int[]
    J = Int[]
    
    for (w1,w2) in zip(W1,W2)
        idx1 = searchsortedfirst(W01,w1)
        idx2 = searchsortedfirst(W02,w2)
        push!(I, idx1)
        push!(J, idx2)
    end
    
    N = length(I)
    
    edist = sparse(I,J,1,length(W01),length(W02),+)
    W01, W02, edist
end

function empdist2Dind(W1::AbstractArray{V1,1}, W2::AbstractArray{V2,1}) where {V1<:Real,V2<:Real}
  
  W01, count1 = countunique(W1)
  W02, count2 = countunique(W2)
  dist = count1 .* count2'
  return W01, W02, dist
end

function concretizearray(inp::Array)
  if 0 == length(inp)
    return inp
  end
  
  newT = typeof(inp[1])
  return map( newT, inp)
end

function trapzdx!(dx::AbstractVector{D}, x::AbstractVector{F}; check_sorted=true) where {F<:Real, D<:Real} 
  @assert length(x) >= 2
  if(check_sorted)
    @assert issorted(x)    
  end
  
  N = length(x)
  resize!(dx, N)

  for i = 2:(N-1)
    dx[i] = (x[i+1] - x[i-1])/2
  end
  dx[1] = (x[2] - x[1])/2
  dx[end] = (x[end] - x[end-1])/2  

  dx
end

function trapzdx(x::AbstractVector{F}; check_sorted=true) where {F<:Real, D<:Real} 
  dx = Vector{F}(undef, length(x))
  trapzdx!(dx, x; check_sorted=check_sorted)
  dx
end

function float_bisection(x1::Float64, x2::Float64)::Float64
    u1 = reinterpret(UInt64, x1)
    u2 = reinterpret(UInt64, x2)
    uM = UInt64(round(u1 //2 + u2//2))
    reinterpret(Float64, uM)
end

function find_crit(x1::Float64, x2::Float64, fun)::Float64
    
    y1 = fun(x1) != 0
    y2 = fun(x2) != 0
    
    for i in 1:64 #max bits to bound   
        if(y1 == y2)
            return NaN
        end
        
        xm = float_bisection(x1,x2)
        
        if(x1==xm) || (x2==xm)
            return xm
        end
        
        ym = fun(xm) != 0
        if (y1==ym)
            x1, y1 = xm, ym
        else
            x2, y2 = xm, ym
        end
        #println("iter=$i, x1=$x1, x2=$x2")
    end
    @warn "no convergence reached"
    return x1
end



end