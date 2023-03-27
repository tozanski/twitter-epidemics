module TwitterCascades

export CASCADE_ALGO
using JLD2
using DataStructures
using Distributions
using Distributed
using SparseArrays
using SharedArrays

@enum CASCADE_ALGO single multitry multi

function cascade_single(
  adj::SparseMatrixCSC{E,V},
  start::Integer, 
  time_distribution::UnivariateDistribution{Support},
  alpha::Real,
  lambda::Real,
  forward::Bool) where {V<:Integer,E,Support}
  
  N, M = size(adj)
  @assert N == M
  @assert start > 0
  @assert start <= N
  @assert alpha <= 1
  @assert alpha >= 0
  @assert lambda >= 0
  
  @assert !forward
  if forward
    adj = copy(adj')
  end
  adj::SparseMatrixCSC{E,V}
  
  tried = falses(N)
  retweeted = falses(N)
    
  is_initiator = true
  queue = [ Pair{Float32,V}(0.0, start) ]
  
  while !isempty(queue)
    recieved_time, vertex = heappop!(queue)
    
    if tried[vertex]
      continue
    end
    tried[vertex] = true
    
    # compute retweeting probability
    retweet_time = recieved_time + rand(time_distribution)
    retweet_prob = alpha * exp(-lambda * retweet_time)
    
    # sample whether vertex retweets (initiator always do)
    is_retweeting = is_initiator || (retweet_prob > rand())
    is_initiator = false
    
    if !is_retweeting
      continue
    end
    
    retweeted[vertex] = true
    
    @inbounds firstedge = adj.colptr[vertex]
    @inbounds lastedge = adj.colptr[vertex+1] - 1
    for edge in firstedge : lastedge
      @inbounds neighbor = adj.rowval[edge]
      
      # if the node has already retweeted we do not allow it to tweet anymore
      @inbounds if tried[neighbor]
        continue
      end
      
      heappush!(queue, Pair{Float32,V}(retweet_time, neighbor))
    end
    
  end
  
  return retweeted
end

function cascade_multitry(
  adj::SparseMatrixCSC{E,V},
  start::Integer, 
  time_distribution::UnivariateDistribution{Support},
  alpha::Real,
  lambda::Real,
  forward::Bool) where {V<:Integer,E,Support}#::BitArray{1}
  
  N, M = size(adj)
  @assert N == M
  @assert start > 0
  @assert start <= N
  @assert alpha <= 1
  @assert alpha >= 0
  @assert lambda >= 0
  
  @assert !forward
  if forward
    adj = copy(adj')
  end
  adj::SparseMatrixCSC{E,V}
  
  retweeted = falses(N)
  is_initiator = true
  queue = [ Pair{Float32,V}(0.0, start) ]
  
  while !isempty(queue)
    recieved_time, vertex = heappop!(queue)
    
    if retweeted[vertex]
      continue
    end
    
    # compute retweeting probability
    retweet_time = recieved_time + rand(time_distribution)
    retweet_prob = alpha * exp(-lambda * retweet_time)
    
    # sample whether vertex retweets (initiator always do)
    is_retweeting = is_initiator || (retweet_prob > rand())
    is_initiator = false
    
    if !is_retweeting
      continue
    end
    
    retweeted[vertex] = true
    
    @inbounds firstedge = adj.colptr[vertex]
    @inbounds lastedge = adj.colptr[vertex+1] - 1
    for edge in firstedge : lastedge
      @inbounds neighbor = adj.rowval[edge]
      
      # if the node has already retweeted we do not allow it to tweet anymore
      @inbounds if retweeted[neighbor]
        continue
      end
      
      # 
      heappush!(queue, Pair{Float32,V}(retweet_time, neighbor))
    end
    
  end
  
  return retweeted
  
end



function cascade(
  algo::CASCADE_ALGO,
  adj::SparseMatrixCSC{E,V},
  start::Integer, 
  time_distribution::UnivariateDistribution{Support},
  alpha::Real,
  lambda::Real,
  forward::Bool) where {V<:Integer,E,Support}
  
  if algo == single
    return cascade_single(adj, start, time_distribution, alpha, lambda, forward)
  elseif algo == multitry
    return cascade_multitry(adj, start, time_distribution, alpha, lambda, forward)
  else
    error("unknown algorithm $algo")
  end

end


function cascadesize(
  algo::CASCADE_ALGO,
  adj::SparseMatrixCSC{E,V}, 
  start::Integer,
  time_distribution::UnivariateDistribution{S1},
  alpha_distribution::UnivariateDistribution{S2},
  lambda_distribution::UnivariateDistribution{S3},
  ;
  forward::Bool=true)::V where {V<:Integer,E,S1,S2,S3}
  
  sum(cascade(algo, adj, start, time_distribution, rand(alpha_distribution), rand(lambda_distribution), forward))
end


function cascadesizes(
  repeats::Integer,
  adj::SparseMatrixCSC{E,V},
  time_distribution::UnivariateDistribution{S1},
  alpha_distribution::UnivariateDistribution{S2},
  lambda_distribution::UnivariateDistribution{S3},
  ;
  forward::Bool=true)::Vector{V} where {V<:Integer,E,S1,S2,S3}

  if forward
    adj = copy(adj')
  end
  adj::SparseMatrixCSC{E,V}
  N::V = size(adj,1)
  
  start_distribution::DiscreteUniform = DiscreteUniform(V(1),N)
  
  local sizes = Vector{V}(repeats)
  local start::V
  local sz::V
  for i in 1:repeats
    start = V(rand( start_distribution ))
    sz = cascadesize(adj,start,time_distribution, alpha_distribution, lambda_distribution, false)
    sizes[i] = V(sz)
  end
  return sizes
end


function cascadesizes_parallel(
  algo::CASCADE_ALGO,
  repeats::Integer,
  adj::SparseMatrixCSC{E,V},
  time_distribution::UnivariateDistribution{S1},
  alpha_distribution::UnivariateDistribution{S2},
  lambda_distribution::UnivariateDistribution{S3},
  ;
  forward::Bool=true) where {V<:Integer,E,S1,S2,S3}

  if forward
    adj = copy(adj')
  end
  adj::SparseMatrixCSC{E,V}
  N::V = size(adj,1)
  
  start_distribution::DiscreteUniform = DiscreteUniform(V(1),N)
  
  sizes = SharedArray{V}(repeats)
  #sizes = Array{V}(undef, repeats)
  
  starts = rand(start_distribution, repeats)
  @sync @distributed for (idx, start) in collect(enumerate(starts))
  #Threads.@threads for (idx, start) in enumerate(starts)
    local sz::V = cascadesize(algo, adj, start, time_distribution, alpha_distribution, lambda_distribution, forward=false)
    sizes[idx] = V(sz)
  end
  
  
  return Array(sizes)
  
end


function tree_cascade_size(
  branching_dist::UnivariateDistribution{Discrete},
  time_dist::UnivariateDistribution, 
  alpha::Real, 
  lambda::Real,
  stack::Vector{ Pair{Float64,Int64} } = Pair{Float64,Int64}[]
)
  

  StackValue = Pair{Float64,Int64}
  
  #stack = StackValue[]
  empty!(stack)
  push!(stack, StackValue(0.0, 1) )
  result = 1

  while !isempty(stack)
    
    received_time, children_count = pop!(stack)
         
    if 0 == children_count
      continue
    end 
    # put the same thing on the stack with children node decreased
    push!(stack, StackValue(received_time, children_count-1) )  
    
    # compute retweeting probability
    retweet_time = received_time + rand(time_dist)
    retweet_prob = alpha * exp(-lambda * retweet_time)
    
    # sample whether vertex retweets (initiator always do)
    is_retweeting = retweet_prob > rand()
    
    if result >= 10^6
      @warn "reached $result retweets, the stack length is $(length(stack))"
      return result
    end    
    
    # if not retweeting close this branch
    if !is_retweeting
      continue
    end
    
   
    
    result += 1
    offspring_count = rand(branching_dist)
    push!(stack, StackValue(retweet_time, offspring_count))
  end
  
  return result
  
end

function tree_cascade_sizes(
  branching_dist::UnivariateDistribution{Discrete},
  time_dist::UnivariateDistribution, 
  alpha_dist::UnivariateDistribution, 
  lambda_dist::UnivariateDistribution,
  N::Integer
)
  
  stack = Pair{Float64,Int64}[]
  
  sizes = Vector{Int}(undef, N)
  
  for i in 1:N
    alpha = rand(alpha_dist)
    lambda = rand(lambda_dist)
    sizes[i] = tree_cascade_size(branching_dist, time_dist, alpha, lambda)
  end
  
  return sizes
  
end

function computecascades!(
  dict, 
  adj::SparseMatrixCSC{E,V},
  time_dist, 
  alpha_dist, 
  niters::Integer, 
  beta_dist = Beta(2,1), 
  lognormal_dist = LogNormal(0,1.5)
) where {V<:Integer, E}
  
  dict["time_dist"] = time_dist
  dict["alpha_dist"] = alpha_dist
  dict["niters"] = niters  
  dict["beta_dist"] = beta_dist
  dict["lognormal_dist"] = lognormal_dist
  @time dict["cascade_single_beta"] = TwitterCascades.cascadesizes_parallel(
    TwitterCascades.single, niters, adj, time_dist, alpha_dist, beta_dist, forward = true);
  
  @time dict["cascade_single_lognormal"] = TwitterCascades.cascadesizes_parallel(
    TwitterCascades.single, niters, adj, time_dist, alpha_dist, lognormal_dist, forward = true);
  
  @time dict["cascade_multitry_beta"] = TwitterCascades.cascadesizes_parallel(
    TwitterCascades.multitry, niters, adj, time_dist, alpha_dist, beta_dist, forward = true);
  
  @time dict["cascade_multitry_lognormal"] = TwitterCascades.cascadesizes_parallel(
    TwitterCascades.multitry, niters, adj, time_dist, alpha_dist, lognormal_dist, forward = true);
end

function computecascades(
  datasetname::AbstractString, 
  filename::AbstractString, 
  adj::SparseMatrixCSC{E,V}, 
  time_dist::UnivariateDistribution, 
  alpha_dist::UnivariateDistribution, 
  niters::Integer, 
  beta_dist::UnivariateDistribution, 
  lognormal_dist::UnivariateDistribution
) where {V<:Integer,E}

  jldopen(filename,"w") do file
    computecascades!(file, adj, time_dist, alpha_dist, niters, beta_dist, lognormal_dist)
    file["datasetname"] = datasetname
  end
  nothing
end

end #module TwitterCascades