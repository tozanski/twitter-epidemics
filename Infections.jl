module Infections

using StatsBase
using BJRG
using SparseArrays


export startsoutbreak
export reachesgiant
export tarjan
export infectionsize
export randinfectionsize

export outbreakprob
export prevalenceprob
export sccprob

function infectionsize(adj::SparseMatrixCSC{E,V}, target::Integer; forward=false, max_discoveries::Integer = V(size(adj,1)), discovery_map::AbstractArray{Bool,1} = zeros(Bool,size(adj,1)), stack::AbstractArray{V,1} = zeros(V,size(adj,1)))::V where {V<:Integer,E}
  N, M = size(adj)
  @assert N == M

  @assert N+1 == length(adj.colptr)
  @assert N <= length(discovery_map)
  @assert N <= length(stack)

  if forward
    adj = copy(adj')
  end
  adj::SparseMatrixCSC{Bool,V}

  stack_ptr = 1                                   # init stack
  stack[1] = target

  discovery_counter = V(1)
  fill!(discovery_map,false)                      # init visit array
  discovery_map[target] = true

  while 0 < stack_ptr # until stack is not empty

    @inbounds v1 = stack[stack_ptr]               # get number of top vertex
    stack_ptr -= 1                                # pop stack

    @inbounds colptr1 = adj.colptr[v1]            # get start of data for top vertex
    @inbounds colptr2 = adj.colptr[v1+1] - 1      # get end of data for top vertex

    @simd for e = colptr1 : colptr2               # for every edge going out from v1
      @inbounds v2 = adj.rowval[e]                # get other node

      @inbounds was_discovered = discovery_map[v2]
      if !was_discovered                             # if v2 was not yet discovered
        @inbounds discovery_map[v2] = true        # mark as discovered
        discovery_counter += 1          # increment counter

        if discovery_counter >= max_discoveries
          return discovery_counter
        end

        stack_ptr += 1
        @inbounds stack[stack_ptr] = v2           # push v2 to the stack
      end
    end
  end
  return discovery_counter
end

function giant_reachability(
  adj::SparseMatrixCSC{Bool,V};
  giant = V(size(adj,1)),
  targets::AbstractArray{V,1} = V(1):V(size(adj,1)),
  as_fraction::Bool = true,
  forward::Bool = false ) where {V<:Integer}

  if forward
    adj = copy(adj')
  end
  adj::SparseMatrixCSC{Bool,V}

  N = V(size(adj,1))

  discovery_map = Vector{Bool}(N)
  stack = Vector{V}(N)

  reaching_giant = V(0)

  if isa( giant, Integer)
    giant_size = V(giant)
    for target in targets
      visited_count = reachability( adj, target, giant_size, discovery_map, stack)
      reaching_giant += visited_count >= giant_size
    end
  elseif isa( giant, AbstractArray) && eltype(giant) <: Integer && ndims(giant) == 1
    @assert length(targets) == length(giant)
    for i in 1:length(targets)
      target = V(targets[i])
      giant_size = V(giant[i])
      visited_count = reachability( adj, target, giant_size, discovery_map, stack)
      reaching_giant += visited_count >= giant_size
    end
  else
    error("giant must be an integer or an array of integers")
  end

  if as_fraction
    return reaching_giant/length(targets)
  end
  return reaching_giant
end


function giant_reachability(adj::SparseMatrixCSC{Bool,V}, transprob::Real;
  targets::AbstractArray{V,1} = V(1):V(size(adj,1)),
  hops::Integer = 3,
  as_fraction::Bool = true,
  forward::Bool = false ) where {V<:Integer}

  fadj = filter_graph(adj,transprob)

  if !forward
    fadj = copy(fadj')
  end
  adj::SparseMatrixCSC{E,V}

  indegree = vec(sum(fadj,1))
  avgdegree = Int(ceil(mean(indegree)))

  giant_sizes = (hops * indegree * avgdegree)[targets] + avg_degree + 1

  giant_reachability(fadj, giant = giant_sizes, targets = targets, as_fraction = as_fraction, forward = false )
end

function tarjan(adj::SparseMatrixCSC{E,V}) where {V<:Integer,E}
  N, M = size(adj)
  @assert N == M


  depthstack = Tuple{V,V}[]

  componentstack = V[]
  incomponent = zeros(Bool,N)

  time::V = V(1)
  times = zeros(V,N)
  low = zeros(V,N)

  componentorder = 1
  componentindices = zeros(V,N)

  for u in V(1):N
    if 0 != times[u]
      continue
    end
    push!(depthstack, (u,V(0)))

    while !isempty(depthstack)
      parent, linkno = pop!(depthstack)
      firstedge = adj.colptr[parent]
      lastedge = adj.colptr[parent+1] - 1
      thisedge = firstedge + linkno
      if 0 == linkno
        times[parent] = time
        low[parent] = time
        push!( componentstack, parent )
        incomponent[parent] = true
        time += 1
      else
        prevchild = adj.rowval[thisedge-1]
        if incomponent[prevchild]
          low[parent] = min(low[parent], low[prevchild])
        end
      end

      if thisedge <= lastedge
        neighbor = adj.rowval[thisedge]
        push!(depthstack, (parent,linkno+1) )
        if 0 == times[neighbor]
          push!( depthstack, (neighbor,V(0)) )
        elseif incomponent[neighbor]
          low[parent] = min(low[parent], times[neighbor])
        end
      else
        if low[parent] == times[parent]
          while last(componentstack) != parent
            w = pop!(componentstack)
            incomponent[w] = false
            componentindices[w] = componentorder
          end
          w = pop!(componentstack)
          incomponent[w] = false
          componentindices[w] = componentorder
          componentorder += 1
        end
      end
    end
  end
  return componentindices
end

function inlargestscc(adj::SparseMatrixCSC{E,V}) where {V<:Integer,E}
    scc = tarjan(adj)
    largest = mode(scc)
    return scc .== largest
end

function inlargestscc(adj::SparseMatrixCSC{E,V}, scc::AbstractArray{Bool,1}) where {V<:Integer,E}
    largest = mode(scc)
    return scc .== largest
end


function startsoutbreak(adj::SparseMatrixCSC{E,V}, scc::AbstractArray{V,1}; forward::Bool=true ) where {V<:Integer,E}
  N, M = size(adj)
  @assert N == M

  if !forward
    adj = copy(transpose(adj))
  end
  adj::SparseMatrixCSC{E,V}

  largest = mode(scc)
  reachable = scc .== largest

  stack = V[]

  for s in V(1):N
    if !reachable[s]
      continue
    end

    firstedge = adj.colptr[s]
    lastedge = adj.colptr[s+1] - 1

    for edge in firstedge : lastedge
      neighbor = adj.rowval[edge]
      if !reachable[neighbor]
        reachable[neighbor] = true
        push!(stack, neighbor)
      end
    end
  end

  while !isempty(stack)
    v = pop!(stack)
    firstedge = adj.colptr[v]
    lastedge = adj.colptr[v+1] - 1

    for edge in firstedge : lastedge
      neighbor = adj.rowval[edge]
      if !reachable[neighbor]
        reachable[neighbor] = true
        push!(stack, neighbor)
      end
    end
  end

  return reachable
end

function startsoutbreak(adj::SparseMatrixCSC{E,V}; forward::Bool=true ) where {V<:Integer,E}
  startsoutbreak(adj, tarjan(adj), forward=forward)
end

function infectionsize(adj::SparseMatrixCSC{E,V}; forward::Bool=true, asfraction::Bool=false) where {V<:Integer,E}
  N, M = size(adj)
  @assert N == M

  scc = tarjan(adj)
  inscc = scc .== mode(scc)
  allwaysinoutbreak = startsoutbreak(adj, scc, forward=!forward)
  minsize = V(sum(allwaysinoutbreak))

  if forward
    adj = copy(adj')
  end
  adj::SparseMatrixCSC{E,V}

  stack = V[]
  sizes = zeros(V,N)
  visited = BitVector(N)

  for s in V(1):N

    if inscc[s]
      sizes[s] = minsize
      continue
    end

    push!(stack, s)
    fill!(visited, false)
    visited[s] = true
    size = 1
    outbreak = false

    while !isempty(stack)
      v = pop!(stack)
      firstedge = adj.colptr[v]
      lastedge = adj.colptr[v+1] - 1

      for edge in firstedge : lastedge
        neighbor = adj.rowval[edge]

        if inscc[neighbor]
          outbreak = true
        elseif !visited[neighbor]
          size += 1
          visited[neighbor] = true
          push!(stack, neighbor)
        end
      end
    end

    sizes[s] = outbreak ? sum(visited | allwaysinoutbreak) : size
  end
  resarr = sizes
  if asfraction
    resarr = sizes .// N
  end
end

function randinfectionsize(adj::SparseMatrixCSC{E,V}, transprob::Real; forward::Bool=true, asfraction::Bool=true) where {V<:Integer,E}
  fadj = filter_graph(adj,transprob)
  return infectionsize(fadj, forward=forward, asfraction=asfraction)
end


function outbreakprob(adj::SparseMatrixCSC{E,V}, transprob::Real) where {V<:Integer,E}
  mean(startsoutbreak(filter_graph(adj,transprob), forward=true) .// V(1))
end

function prevalenceprob(adj::SparseMatrixCSC{E,V}, transprob::Real) where {V<:Integer,E}
  mean(startsoutbreak(filter_graph(adj,transprob), forward=false) .// V(1))
end

function sccprob(adj::SparseMatrixCSC{E,V}, transprob) where {V<:Integer,E}
  sccno = tarjan(filter_graph(adj,transprob))
  largest = mode(sccno)
  mean( (sccno .== largest) .// V(1) )
end

function resistant(adj::SparseMatrixCSC{E,V}, infected::AbstractArray{Bool,1}; forward=true) where {V<:Integer,E}

  N = size(adj,1)
  M = size(adj,2)
  @assert N==M
  @assert N == length(infected)

  if forward
    adj = copy(adj')
  end
  adj::SparseMatrixCSC{E,V}

  resistant = falses(N)

  for (v,is_infected) in enumerate(infected)
    if !is_infected
      continue
    end

    firstedge = adj.colptr[v]
    lastedge = adj.colptr[v+1] - 1

    for edge in firstedge : lastedge
      neighbor = adj.rowval[edge]
      if !infected[neighbor]
        resistant[neighbor] = true
      end
    end
  end


  return resistant
end


end
