module Twitter

using BJRG
using Distributed
using JLD2
using Infections
using TomekUtils
using StatsBase
using SparseArrays
using Random
using LinearAlgebra

export follower_estimate
export follower_estimate_scalar
export follower_similar
export follower_similar_scalar
export bjrsolve
export bjrsolve_ind
#export bjrsolve_undirected

export computeinfections!
export computeinfections

function follower_estimate(adj::AbstractArray{T,2}) where T<:Real
    @assert size(adj,1) == size(adj,2)
    d1 = vec(sum(adj,dims=2))
    d2 = vec(sum(adj,dims=1))
    m = mean(d1)
    w1 = d1/sqrt(m)
    w2 = d2/sqrt(m)
    return w1,w2,1.0
end

function follower_estimate_scalar(adj::AbstractArray{T,2},gamma::Real) where {T<:Real}
    @assert size(adj,1) == size(adj,2)
    local d1 = vec(sum(adj,dims=2))
    local d2 = vec(sum(adj,dims=1))

    local f1 = mean(d1.^gamma).^(-1/(1+gamma))
    local f2 = mean(d2.^(1/gamma)).^((-gamma)/(gamma+1))

    local w1 = d1 * f1
    local w2 = d2 * f2

    local w = (w1 + w2.^(1/gamma))/2
    return w
end

function likelihood_scalar(adj::SparseMatrixCSC{E,V},w::AbstractArray{W,1},gamma::Float64) where {V<:Real,E,W}
    local N = length(w)
    local wg = w.^gamma
    local le = log.( nonzeros( Diagonal(w) * adj * Diagonal(wg) ) )
    local lb = sum(w) * sum(wg)/N
    return sum(le) - lb - log(N)*nnz(adj)
end

function follower_estimate_scalar(adj::AbstractArray{T,2},gamma_min::Real, gamma_max::Real; reltol::Real=0, maxiters::Integer=100) where {T<:Real}

  local g1 = log(gamma_min)
  local g4 = log(gamma_max)
  local g2 = g4 - (g4 - g1)/MathConstants.golden
  local g3 = g1 + (g4 - g1)/MathConstants.golden

  @assert 0 <= exp(g1)
  @assert g1 <= g2
  @assert g2 <= g3
  @assert g3 <= g4

  local w1 = follower_estimate_scalar(adj, exp(g1))
  local l1 = likelihood_scalar(adj, w1, exp(g1))

  local w2 = follower_estimate_scalar(adj, exp(g2))
  local l2 = likelihood_scalar(adj, w2, exp(g2))

  local w3 = follower_estimate_scalar(adj, exp(g3))
  local l3 = likelihood_scalar(adj, w3, exp(g3))

  local w4 = follower_estimate_scalar(adj, exp(g4))
  local l4 = likelihood_scalar(adj, w4, exp(g4))

  local converged::Bool
  converged = false

  for iter in 1:maxiters
    #println((g1, g2, g3, g4))
    #println((exp(g1), exp(g2), exp(g3), exp(g4)))
    #println((l1, l2, l3, l4))

    if reltol==0 && nextfloat(exp(g2))>=exp(g3) || reltol!=0 && g3-g2 < reltol
      converged = true
      break

    end

    if l2 > l3  # maximum between l1 and l3, g1 stays, g4 can be discarded
      g4 = g3; w4 = w3; l4 = l3
      g3 = g2; w3 = w2; l3 = l2
      g2 = g4 - (g4-g1)/MathConstants.golden
      w2 = follower_estimate_scalar(adj, exp(g2))
      l2 = likelihood_scalar(adj, w2, exp(g2))
    else # maximum between l2 and l4, g1 can be discarded, g4 stays
      g1 = g2; w1 = w2; l1 = l2
      g2 = g3; w2 = w3; l2 = l3
      g3 = g1 + (g4-g1)/MathConstants.golden
      w3 = follower_estimate_scalar(adj, exp(g3))
      l3 = likelihood_scalar(adj, w3, exp(g3))
    end
  end
  if !converged
    @warn "estimation did not converge in $maxiters iterations"
  end
  return l2 > l3 ? (w2, exp(g2)) : (w3, exp(g3))
end

function follower_kernel(x,y)
    @inbounds v = x[1]*y[2]
    return v
end

function follower_similar(adj::AbstractArray{T,2}; shuffled::Bool = false) where {T<:Real}
    w1,w2,m = follower_estimate(adj)

    local m::Float64
    w1 *= m

    if shuffled
      shuffle!(w2)
    end
    weights = collect(zip(w1,w2))

    adj2 = generate_bjr_graph(weights, follower_kernel);

    return adj2
end

function follower_similar_scalar(adj::SparseMatrixCSC{E,V}, gamma_min::G, gamma_max::G) where {V<:Integer,E,G<:Real}
  #const gamma::Float64
  #const w::Array{Float64,1}

  w0, gamma0 = follower_estimate_scalar(adj, gamma_min, gamma_max)
  gamma::Float64 = gamma0
  w::Vector{Float64} = w0

  function kernel(x,y)::Float64
    x*y.^gamma
  end
  #const kernel = (x,y) -> x * y^gamma
  adj2 = generate_bjr_graph(w, kernel)
  return adj2
end

function bjrsolve(w1::AbstractVector{W1}, w2::AbstractVector{W2}, edist::AbstractArray{C,2}, c::Real) where {W1<:Real,W2<:Real,C<:Real}
  MAX_ITERS = 100
  # double buffering scheme
  rho = ones(W1, size(w1))
  old_rho = ones(W1, size(w1))

  v = Array{W1,1}(undef, length(w2))

  m = c / sum(edist)

  for i in 1:MAX_ITERS
    mul!(v, transpose(edist), rho)
    #v = edist' * rho
    f = dot(v, w2) * m

    # swap buffers
    rho, old_rho = old_rho, rho

    rho .= 1 .- exp.(-f .* w1)
#    for i in 1:length(rho)
#      rho[i] = 1 - exp(-f*w1[i])
#    end

    if rho == old_rho
      return rho
    end

  end
  @warn "no convergence reached after $MAX_ITERS iterations"
  return rho
end

function bjrsolve(w1::AbstractVector{W1}, w2::AbstractVector{W2}, dist1::AbstractVector{C1}, dist2::AbstractVector{C2}, c::Real) where {W1<:Real,W2<:Real,C1<:Integer,C2<:Integer}

  MAX_ITERS = 100
  @assert c >= 0

  rho = ones(W1, size(w1) )
  old_rho = ones(W1, size(w1) )
  f = c / sum(dist1) / sum(dist2)

  for i in 1:MAX_ITERS
    rho, rho_old = old_rho, rho
    v = dot(rho,dist1) * dot(w2,dist2) * f
    rho .= 1 - exp.(-v * w1 )
    if rho_old == rho
      return rho
    end
  end
  @warn "no convergence reached in $MAX_ITERS iterations"
  return rho
end

function bjrsolveind(w1::AbstractVector{W1}, w2::AbstractVector{W2}, dist1::AbstractVector{C1}, dist2::AbstractVector{C2}, c::Real) where {W1<:Real,W2<:Real,C1<:Integer,C2<:Integer}

  MAX_ITERS = 100
  @assert c >= 0

  rho = ones(W1, size(w1) )
  old_rho = ones(W1, size(w1) )
  f = c / sum(dist1) / sum(dist2)

  for i in 1:MAX_ITERS
    rho, rho_old = old_rho, rho
    v = dot(rho,dist1) * dot(w2,dist2) * f
    rho = 1 - exp(-v * w1 )
    if rho_old == rho
      return rho
    end
  end
  @warn "no convergence reached in $MAX_ITERS iterations"
  return rho
end

function bjrresistant(rho::AbstractVector{R}, w1::AbstractVector{W1}, w2::AbstractVector{W2}, dist::AbstractArray{C,2}, p::Real) where {R<:Real,W1<:Real,W2<:Real,C<:Real}

  N,M = size(dist)
  @assert length(w1) == N
  @assert length(w2) == M
  @assert length(rho) == M

  dist1 = vec(sum(dist,dims=2)) / sum(dist)
  dist2 = vec(sum(dist,dims=1)) / sum(dist)

  @assert length(dist1) == N
  @assert length(dist2) == M

  lambda0 = sum( rho .* dist2 .* w2)
  lambda = w1 .* lambda0 # *p ?


  @assert length(lambda) == N

  prob = 1 .- exp.(-lambda)

  #(1-rho) #* dist / sum(dist)
  dot(prob, dist * (1 .- rho)) / sum(dist)



#    @assert length(w1) == size(dist,1)
#    @assert length(w2) == size(dist,2)
#    @assert length(rho) == length(w1)
#
#    distnorm = sum(dist)
#
#    v0 = 1./dot( vec( sum(dist,2)), w1 )
#    v1 = dot( (1-rho), vec(sum(dist,2)) )
#    v2 = dot( vec(sum(dist,2)), rho .* w1) / distnorm

#    return (1-p) .* v0 .* v1 .* v2
    #sum( (1-p) * sum( (1-rho) .* rho' ) ) * (1-p)/N
end

function iterative_step_ind(rho::AbstractArray{V,1}, w1::AbstractArray{V,1}, w2::AbstractArray{V,1}, c::Real, dist1 = ones(w1), dist2 = ones(w2)) where {V<:Real}
  v = dot(rho,dist1) * dot(w2,dist2) / sum(dist1) / sum(dist2)
  1 - exp.(-c * v * w1 )
end

function iterative_step_ind(rho::AbstractArray{V,1}, w1::AbstractArray{V,1}, w2::AbstractArray{V,1}, c::Real, dist1::Nothing, dist2::Nothing) where {V<:Real}

  v = mean(rho) * mean(w2)
  1 - exp.(-c * v * w1 )
end

function bjr_solve_ind2(w1::AbstractArray{V,1}, w2::AbstractArray{V,1}, c::Real, dist1 = nothing, dist2 = nothing) where {V<:Real}

  rho = ones(w1)
  for i in 1:100
    rho_old = rho
    rho = iterative_step_ind(rho, w1, w2, c, dist1, dist2)
    if rho_old == rho
      return rho
    end
  end
  return rho
end



function simulateinfections(adj::SparseMatrixCSC{E,V},p::Real) where {V<:Real,E}
  filtered = BJRG.filter_graph(adj,p)
  scc = Infections.tarjan(filtered)
  N = V(length(scc))

  core = sum(scc .== StatsBase.mode(scc) ) // N
  initiators = Infections.startsoutbreak(filtered, scc, forward=true)
  outbreak = sum( initiators ) // N
  infected =  Infections.startsoutbreak(filtered, scc, forward=false)
  prevalence = sum(infected) // N
  resistant = sum( Infections.resistant(adj, infected, forward=false) ) // N

  return core, outbreak, prevalence, resistant
end

function simulateinfections(adj::SparseMatrixCSC{E,V}, transprob::AbstractArray{P,1}) where {V<:Real,P<:Real,E}

  N = length(transprob)

  results = @distributed hcat for p in transprob
    core, outbreak, prevalence, resistant = simulateinfections(adj, p)
    vcat(core, outbreak, prevalence, resistant)
  end

  cores = results[1,:]
  outbreaks = results[2,:]
  prevalences = results[3,:]
  resistances = results[4,:]
  contacts = prevalences + resistances

  return cores, outbreaks, prevalences, resistances, contacts
end

function bjrinfections(w1::AbstractArray{W1,1}, w2::AbstractArray{W2,1}, edist::AbstractArray{C,2},p::Real) where {W1<:Real, W2<:Real, C<:Real}
  edistnorm = sum(edist)

  outbreak_rho = bjrsolve(w1,w2,edist,p)
  outbreak_cnt = vec( sum(edist,dims=2) ) / edistnorm
  outbreak_tot = dot( outbreak_rho, outbreak_cnt )

  prevalence_rho = bjrsolve(w2,w1,edist', p)
  prevalence_cnt = vec( sum(edist,dims=1) ) / edistnorm
  prevalence_tot = dot( prevalence_rho, prevalence_cnt )

  resistance_tot = bjrresistant(prevalence_rho, w1, w2, edist, p)

  core_tot = dot(outbreak_rho, edist * prevalence_rho) / edistnorm
  return core_tot, outbreak_tot, prevalence_tot, resistance_tot
end

function bjrinfections(w1::AbstractArray{W1,1}, w2::AbstractArray{W2,1}, edist::AbstractArray{C,2},transprob::AbstractVector{P}) where {W1<:Real, W2<:Real, C<:Real, P<:Real}

  outbreaks = Float64[]
  prevalences = Float64[]
  cores = Float64[]
  resistances = Float64[]

  for p in transprob
    core, outbreak, prevalence, resistance = bjrinfections(w1, w2, edist, p)
    push!(cores, core)
    push!(prevalences, prevalence)
    push!(outbreaks, outbreak)
    push!(resistances, resistance)
  end

  contacts = resistances + prevalences

  return cores, outbreaks, prevalences, resistances, contacts
end

function computeinfections!(dict, transprob::AbstractArray{P,1}, adj, adj_est, adj_shu, adj_sca, gamma) where P
  dict["transprob"] = transprob
  dict["gamma"] = gamma

  w1,w2,m = follower_estimate(adj)
  w0 = follower_estimate_scalar(adj, gamma)

  w01, w02, dist = TomekUtils.empdist2D(w1,w2*m)
  @time dict["scc_theo"], dict["outbreak_theo"], dict["prevalence_theo"], dict["resistance_theo"], dict["contact_theo"] = bjrinfections(w01, w02, dist, transprob)

  w01, w02, dist = TomekUtils.empdist2Dind(w1,w2*m)
  @time dict["scc_theo_ind"], dict["outbreak_theo_ind"], dict["prevalence_theo_ind"], dict["resistance_theo_ind"], dict["contact_theo_ind"] = bjrinfections(w01, w02, dist, transprob)

  w01, w02, dist = TomekUtils.empdist2D(w0, w0.^gamma)
  @time dict["scc_theo_sca"], dict["outbreak_theo_sca"], dict["prevalence_theo_sca"], dict["resistance_theo_sca"], dict["contact_theo_sca"] = bjrinfections(w01, w02, dist, transprob)


  @time dict["scc"], dict["outbreak"], dict["prevalence"], dict["resistance"], dict["contact"] = simulateinfections(adj, transprob)
  @time dict["scc_est"], dict["outbreak_est"], dict["prevalence_est"], dict["resistance_est"], dict["contact_est"] = simulateinfections(adj_est, transprob )
  @time dict["scc_shu"], dict["outbreak_shu"], dict["prevalence_shu"], dict["resistance_shu"], dict["contact_shu"] = simulateinfections(adj_shu, transprob )
  @time dict["scc_sca"], dict["outbreak_sca"], dict["prevalence_sca"], dict["resistance_sca"], dict["contact_sca"] = simulateinfections(adj_sca, transprob )


end

function computeinfections(datasetname::AbstractString, fname::AbstractString, transprob, adj, adj_est, adj_shu, adj_sca, gamma)
  jldopen(fname,"w") do file
    computeinfections!(file, transprob, adj, adj_est, adj_shu, adj_sca, gamma)
    file["datasetname"] = datasetname
  end
  nothing
end

end
