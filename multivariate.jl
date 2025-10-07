struct MultivariateInducedDistribution <: ContinuousMultivariateDistribution
    index_set::Vector{<:Tuple}
    N::Int
    polys::Tuple
    PolyObjs::Vector{AbstractOrthoPoly}
    x_Qs::Vector{Vector{Float64}}
    w_Qs::Vector{Vector{Float64}}
    normssqs::Vector{Vector{Float64}}
    Pns::Vector{Matrix{Float64}}
    Vs::Vector{Matrix{Float64}}
    Ws::Vector{Matrix{Float64}}
    Fs::Vector{Matrix{Float64}}
end

Base.length(d::MultivariateInducedDistribution) = length(d.polys)

function Distributions._rand!(::AbstractRNG, d::MultivariateInducedDistribution, res::AbstractVector{<:Real})
    N = length(d.index_set)
    D = length(d)
    n = rand(1:N)
    idx = d.index_set[n]
    for dim in 1:D
        u = rand()
        nidx = idx[dim]
        i = searchsortedfirst(view(d.Fs[dim], nidx+1, :), u) # ∈ {1,...,Q}
        res[dim] = d.x_Qs[dim][i]
    end
    return res
end

Distributions.sampler(d::MultivariateInducedDistribution) = d
Base.eltype(d::MultivariateInducedDistribution) = Float64

function Distributions._logpdf(d::MultivariateInducedDistribution,x::AbstractArray)
    # = dρ(x) = 1/N ∑ pn(x)^2 dμ(x)
    N = length(d.index_set)
    D = length(d)
    res = 0.0
    Pn²s = [evaluate(0:size(d.Pns[j],1)-1, x[j], d.PolyObjs[j]) .^ 2 for j in 1:D]
    for idx in d.index_set
        cur = 1.0
        for j in 1:D
            cur *= Pn²s[j][idx[j]+1] / d.normssqs[j][idx[j]+1]
        end
        res += cur
    end
    for j in 1:D 
        res *= pdf(_dists[d.polys[j]], x[j])
    end
    res /= N
    return log(res)
end

function MultivariateInducedDistribution(polys=(:legendre,:legendre), N=5; index_set=multi_index_set(length(polys), N), Q=500)
    D = length(polys)
    
    all_measures = [_measures[poly]() for poly in polys] 
    PolyObjs = [OrthoPoly(measure, N, Nrec=Q+1) for measure in all_measures]

    x_Qs = [PolyObjs[j].quad.nodes for j in 1:D]
    w_Qs = [PolyObjs[j].quad.weights for j in 1:D]

    x_Q = Vector{Float64}[]
    w_Q = Float64[]
    Q_index_set = multi_index_set(D, Q-1, :tp)
    # for idx in Q_index_set
    #     x = [x_Qs[j][idx[j]+1] for j in 1:D]
    #     w = prod(w_Qs[j][idx[j]+1] for j in 1:D)
    #     push!(x_Q, x)
    #     push!(w_Q, w)
    # end

    # X_Q = reduce(hcat, x_Q)

    # Rename vars
    #N = length(index_set)
    # Q = length(x_Q)

    # Evaluate polynomials and normalize
    Pns = [
        collect(transpose(evaluate(0:N, PolyObjs[d].quad.nodes, PolyObjs[d])))
            for d in 1:D
    ]
    normssqs = [
        (Pns[d] .^ 2) * w_Qs[d] for d in 1:D
    ]
    # Compute V, W, F matrices
    Vs = [
        Pns[d] ./ sqrt.(normssqs[d]) for d in 1:D        
    ]
    Ws = [
        Vs[d] .* Vs[d] .* w_Qs[d]' for d in 1:D
    ]
    for d in 1:D
        Ws[d] ./= sum(Ws[d], dims=2) # For stability, not necessary
    end
    Fs = [
        cumsum(Ws[d], dims=2) for d in 1:D
    ]
    return MultivariateInducedDistribution(index_set, N, Tuple(polys), PolyObjs, x_Qs, w_Qs, normssqs, Pns, Vs, Ws, Fs)
end
