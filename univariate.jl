struct InducedDistribution <: ContinuousUnivariateDistribution
    N::Int
    poly::Symbol
    PolyObj::AbstractOrthoPoly
    x_Q::Vector{Float64}
    w_Q::Vector{Float64}
    normssq::Vector{Float64}
    Pn::Matrix{Float64}
    V::Matrix{Float64}
    W::Matrix{Float64}
    F::Matrix{Float64}
end

function Base.rand(::AbstractRNG, d::InducedDistribution)
    N = d.N
    n = rand(1:(N+1))
    u = rand()
    i = searchsortedfirst(view(d.F, n, :), u) # ∈ {1,...,Q}
    return d.x_Q[i]
end

Distributions.sampler(d::InducedDistribution) = d

function Distributions.logpdf(d::InducedDistribution,x::Real)
    # = dρ(x) = 1/N ∑ pn(x)^2 dμ(x)
    res = 0.0
    Pn² = evaluate(0:d.N, x, d.PolyObj) .^ 2
    for i in 1:(d.N+1)
        res += Pn²[i] / d.normssq[i]
    end
    res *= pdf(_dists[d.poly], x) / (d.N + 1)
    return log(res)
end

function Distributions.cdf(d::InducedDistribution,x::Real)
    f = sum(d.F, dims=1)[1,:] ./ size(d.F, 1)
    i = searchsortedfirst(d.x_Q, x)
    i = min(i, length(f))
    return f[i]
end
function Distributions.quantile(d::InducedDistribution,q::Real)
    # Inverse CDF function
    cdfs = cdf.(d, d.x_Q)
    i = searchsortedfirst(cdfs, q)
    i = min(i, length(cdfs))
    return d.x_Q[i]
end

Base.minimum(d::InducedDistribution) = minimum(d.x_Q)
Base.maximum(d::InducedDistribution) = maximum(d.x_Q)
Distributions.insupport(d::InducedDistribution,x::Real) = minimum(d) <= x <= maximum(x)#insupport(dists[d.poly],x)

function InducedDistribution(poly=:legendre, N=5; Q=500)
    measure = _measures[poly]()
    PolyObj = OrthoPoly(measure, N, Nrec=Q+1)

    x_Q = PolyObj.quad.nodes
    w_Q = PolyObj.quad.weights

    # Evaluate polynomials and normalize
    Pn = collect(transpose(evaluate(0:N, x_Q, PolyObj)))
    # Compute norms of each pn
    normssq = (Pn .^ 2) * w_Q
    # Compute V, W, F matrices
    V = Pn ./ sqrt.(normssq)
    W = V .* V .* w_Q'
    W ./= sum(W, dims=2) # For stability, not necessary
    F = cumsum(W, dims=2)

    return InducedDistribution(N, poly, PolyObj, x_Q, w_Q, normssq, Pn, V, W, F)
end