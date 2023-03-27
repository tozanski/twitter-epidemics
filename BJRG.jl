module BJRG

using Distributions
using Random
using SparseArrays

export generate_bjr_graph
export generate_er_graph
export filter_graph

function generate_bjr_graph(weights::AbstractArray{T,1}, kernel::Function; seed::Integer=0, directed::Bool = true, selfedges::Bool = false) where T

  rng = MersenneTwister(seed)
  N = length(weights)
  I = Int[]
  J = Int[]

  for j in 1:N
    for i in 1:N
      if selfedges && i!=j
        continue
      end

      if !directed && i<j
        continue
      end

      val::Real = kernel(weights[i], weights[j]) / N
      if rand(rng) < val
        push!(I,i)
        push!(J,j)
      end
    end
  end

  return sparse(I,J,true,N,N)
end

function generate_er_graph(N::V, meandeg::Real, rng::AbstractRNG=MersenneTwister(0), selfedges=false) where V <: Integer

  edgenum_dist = Poisson(meandeg*N)
  numel_est = cquantile(edgenum_dist, 10.0^-6) # chance of relocation 1 in million

  colptr = V[]
  sizehint!(colptr,N+1)

  rowval = V[]
  sizehint!(rowval,numel_est)

  skip_dist = Geometric(meandeg/N) # in case of non selfedges we will have slightly more samples they will be compensated later

  push!(colptr,1)
  ptr = 1
  row = 1
  col = 1
  while col <= N

    skip = rand(rng, skip_dist) # how many non-edges are skipped before an edge is inserted

    row += skip + 1

    while row > N
      push!(colptr, ptr)
      row -= N
      col += 1
    end

    if selfedges || row != col # if self edges are not allowed do not place them it is compensated by sampling slightly more edges than required
      push!(rowval, row)
      ptr += 1
    end


  end

  nzval = ones(Bool, length(rowval))
  return SparseMatrixCSC{Bool,V}(N,N,colptr,rowval,nzval)

end

function filter_graph(adj::SparseMatrixCSC{Bool,V}, p::Real; seed = 0)::SparseMatrixCSC{Bool,V} where V <: Integer
  N, M = size(adj)
  @assert N == M
  @assert p <= 1.0
  @assert 0 <= p

  rng = MersenneTwister(seed)

  new_colptr = Vector{V}(undef, N+1)
  new_rowval = V[]
  expected_edges = nnz(adj) * p
  sizehint!(new_rowval, Int(ceil(expected_edges + sqrt(expected_edges)) ))
  new_edgeno = 1

  new_colptr[1] = 1

  for v in V(1):V(N)
    firstedge = adj.colptr[v]
    lastedge = adj.colptr[v+1] - 1

    for edge in firstedge : lastedge
      neighbor = adj.rowval[edge]

      if rand(rng) < p
        push!(new_rowval, neighbor)
        new_edgeno += 1
      end
    end
    new_colptr[v+1] = new_edgeno

  end
  new_nzval = ones(Bool,length(new_rowval))

  return SparseMatrixCSC(N, N, new_colptr, new_rowval, new_nzval)
end

end