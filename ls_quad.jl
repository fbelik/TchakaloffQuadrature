mutable struct LS_Quad{T}
    Mtot::Int
    const d::Union{InducedDistribution,MultivariateInducedDistribution}
    V::Matrix{T}
    Vsq::Matrix{T}
    const wsq::Vector{T}
    const τ²::Vector{T}
    const Gm::Matrix{T}
    Gminv::Union{Matrix{T},Nothing}
    const η::Vector{T}
    const nodes::Vector
    const weights::Vector{T}
    const induced::Bool
end

# Assuming G=I
function least_squares_quad(d::InducedDistribution, M=100; induced=true)
    if induced
        nodes = rand(d, M)
    else
        d0 = _dists[d.poly]
        nodes = rand(d0, M)
    end
    V = evaluate(0:d.N, nodes, d.PolyObj) ./ sqrt.(d.normssq)'
    N = size(V, 2)
    τ² = sum(V .^ 2, dims=2)[:,1] ./ N
    if !induced
        τ² .= 1.0
    end
    Vsq = Matrix{Float64}(undef, M, Int(N*(N+1)/2))
    idx = 1
    for i in 1:N
        for j in i:N
            for m in 1:M
                Vsq[m,idx] = V[m,i] * V[m,j] / τ²[m]
            end
            idx += 1
        end
    end
    wsq = ones(M) ./ M
    Gm = (V' * Diagonal(1 ./ τ²) * V) ./ M
    η = d.V * d.w_Q
    w = zeros(M)
    Gminv = nothing
    try
        w .= (Diagonal(1 ./ τ²) * (V * (Gm \ η))) ./ M
        Gminv = inv(Gm)
    catch SingularException
        w .= (Diagonal(1 ./ τ²) * (V * (pinv(Gm) \ η))) ./ M
        Gminv = nothing
    end
    return LS_Quad(M, d, V, Vsq, wsq, τ², Gm, Gminv, η, nodes, w, induced)
end

function least_squares_quad(d::MultivariateInducedDistribution, M=100; induced=true)
    D = length(d.polys)
    N = length(d.index_set)
    if induced
        nodes = rand(d, M)
    else
        d0s = [_dists[poly] for poly in d.polys]
        nodes = reduce(hcat, rand(d0, M) for d0 in d0s)'
    end
    Vs = [
        evaluate(0:d.N, view(nodes,j,:), d.PolyObjs[j]) ./ sqrt.(d.normssqs[j])'
            for j in 1:D
    ]
    V = ones(M, N)
    for m in 1:M
        for n in 1:N
            idx = d.index_set[n]
            for (ii,i) in enumerate(idx)
                V[m,n] *= Vs[ii][m,i+1]
            end
        end
    end
    τ² = sum(V .^ 2, dims=2)[:,1] ./ N
    if !induced
        τ² .= 1.0
    end  
    Vsq = Matrix{Float64}(undef, M, Int(N*(N+1)/2))
    idx = 1
    for i in 1:N
        for j in i:N
            for m in 1:M
                Vsq[m,idx] = V[m,i] * V[m,j] / τ²[m]
            end
            idx += 1
        end
    end
    wsq = ones(M) ./ M
    η = zeros(N); η[1] = 1;#d.V * d.w_Q
    nodes = [collect(x) for x in eachcol(nodes)]
    Gm = (V' * Diagonal(1 ./ τ²) * V) ./ (M)
    G = I
    w = (Diagonal(1 ./ τ²) * (V * (Gm \ η)))./ M
    Gminv = nothing
    try
        Gminv = inv(Gm)
    catch SingularException
    end
    return LS_Quad(M, d, V, Vsq, wsq, τ², Gm, Gminv, η, nodes, w, induced)
end

function add_node!(q::LS_Quad, x=nothing; nodeidx=nothing, gminvtol=1e-2)
    M = q.Mtot
    N = size(q.V, 2)
    # Sample new node
    if isnothing(x)
        if q.induced
            x = rand(q.d)
        elseif isa(q.d, InducedDistribution)
            x = rand(_dists[q.d.poly])
        else
            x = [rand(_dists[p]) for p in q.d.polys]
        end
    end
    # Compute new row of v
    if isa(q.d, InducedDistribution)
        v = evaluate(0:d.N, x, d.PolyObj) ./ sqrt.(d.normssq)'
        v = view(v, 1, :)
        τ²new = sum(v .^ 2) ./ N
    else
        D = length(d.polys)
        vs = [evaluate(0:d.N, x[j], d.PolyObjs[j]) ./ sqrt.(d.normssqs[j])' 
                for j in 1:D]
        v = ones(N)
        for (i,idx) in enumerate(d.index_set)
            for (ct,j) in enumerate(idx)
                v[i] *= vs[ct][j+1]
            end
        end
        τ²new = sum(v .^ 2) ./ N
    end
    if !q.induced
        τ²new = 1.0
    end
    if isnothing(nodeidx)
        # Add new row to V
        q.V = vcat(q.V, v')
        # Update τ²
        push!(q.τ², τ²new)
        # Add new row to Vsq
        q.Vsq = vcat(q.Vsq, Vector{Float64}(undef, Int(N*(N+1)/2))')
        idx = 1
        for i in 1:N
            for j in i:N
                q.Vsq[end,idx] = v[i] * v[j] / τ²new
                idx += 1
            end
        end
        # Add new element to wsq
        q.wsq .*= M / (M+1)
        push!(q.wsq, 1 / (M+1))
        # Add node
        push!(q.nodes, x)
        # Update weights
        push!(q.weights, 0)
        q.weights .= (Diagonal(1 ./ q.τ²) * (q.V * (q.Gminv * q.η))) ./ (M+1)
    else
        # Insert new row to V
        q.V[nodeidx,:] .= v
        # Update τ²
        q.τ²[nodeidx] = τ²new
        # Insert new row to Vsq
        idx = 1
        for i in 1:N
            for j in i:N
                q.Vsq[nodeidx,idx] = v[i] * v[j] / τ²new
                idx += 1
            end
        end
        # Update wsq
        q.wsq .*= M / (M+1)
        q.wsq[nodeidx] = 1 / (M+1)
        # Insert node
        q.nodes[nodeidx] = x
    end
    # Update Gm
    q.Gm .= M / (M+1) .* q.Gm .+ (v * v') ./ (τ²new * (M+1))
    # Update Gm⁻¹ by Sherman-Morrison
    if !isnothing(q.Gminv)
        prod = q.Gminv * v
        q.Gminv .= (M+1)/M .* (q.Gminv .- (prod * prod') ./ (M*τ²new + dot(v,q.Gminv,v)))
    else
        try
            q.Gminv = inv(q.Gm)
        catch SingularException
            q.Gminv = nothing
        end
    end
    # Update weights
    q.weights .= (Diagonal(1 ./ q.τ²) * (q.V * (q.Gminv * q.η))) ./ (M+1)
    # Check Gminv
    if !isnothing(q.Gminv) && norm(q.Gm * q.Gminv - I) > gminvtol
        try
            q.Gminv .= inv(q.Gm)
        catch SingularException
            q.Gminv = nothing
        end
    end
    q.Mtot += 1
    return
end

function remove_add_node!(q::LS_Quad, x=nothing; gminvtol=1e-2)
    M = q.Mtot
    N = size(q.V, 2)
    # Select node to remove from q
    Q, _ = qr(q.Vsq)
    n = Q[:,end]
    remval, remidx = findmin(abs.(q.wsq ./ n))
    α = q.wsq[remidx] / n[remidx]
    q.wsq .-= α .* n
    # Remove node remidx from q
    new_indices = [1:remidx-1 ; remidx+1:size(q.V,1)]
    # Update weights, and V or τ²
    #q.V .= Diagonal(q.wsq ./ (q.wsq .+ α .* n)) * q.V
    q.τ² ./= (q.wsq ./ (q.wsq .+ α .* n))
    q.weights[new_indices] .= (Diagonal(1 ./ q.τ²[new_indices]) * (view(q.V,new_indices,:) * (q.Gminv * q.η))) ./ M
    # Add new node at remidx
    add_node!(q, x, nodeidx=remidx, gminvtol=gminvtol)
end

function combine_weights(q::LS_Quad{T}) where T
    dict = Dict{eltype(q.nodes),T}()
    inds = Int[]
    for (i,(x,w)) in enumerate(zip(q.nodes,q.weights))
        if x in keys(dict)
            dict[x] += w
        else
            dict[x] = w
            push!(inds, i)
        end
    end
    nodes = collect(keys(dict))
    weights = [dict[x] for x in nodes]
    return LS_Quad{T}(q.d, q.V[inds,:], q.Gm, q.nodes[inds], weights)
end

function PolyChaos.evaluate(q::LS_Quad, f)
    res = 0.0
    for (x,w) in zip(q.nodes,q.weights)
        res += w * f(x)
    end
    return res
end

function concentration(q::LS_Quad)
    return opnorm(q.Gm - I)
end

function Plots.plot(q::LS_Quad; kwargs...)
    if isa(q.d, InducedDistribution)
        scatter(q.nodes, q.weights; kwargs...)
    else
        scatter([x[1] for x in q.nodes], 
                [x[2] for x in q.nodes],
                q.weights; kwargs...)
    end
end

function Plots.plot!(q::LS_Quad; kwargs...)
    if isa(q.d, InducedDistribution)
        scatter!(q.nodes, q.weights; kwargs...)
    else
        scatter!([x[1] for x in q.nodes], 
                [x[2] for x in q.nodes],
                q.weights; kwargs...)
    end
end

function is_positive(q::LS_Quad; tol=0.0)
    mval = -1 * abs(tol)
    for w in q.weights
        if w < mval
            return false
        end
    end
    return true
end

function trial_comparison(d::Union{InducedDistribution,MultivariateInducedDistribution}, M=100, trials=1000; tol=1e-4, doplot=true, pbar=true)
    succ = [0,0]
    for _ in (pbar ? ProgressBar(1:trials) : 1:trials)
        q = least_squares_quad(d, M);
        q0 = least_squares_quad(d, M, induced=false);
        is_positive(q, tol=tol) ? succ[1] += 1 : nothing
        is_positive(q0, tol=tol) ? succ[2] += 1 : nothing
    end
    if doplot
        plt = bar(100 .* succ ./ trials, label=false)
        plot!(plt, xticks=(1:2,("Induced","Standard")), ylabel="Success %")
        plot!(ylims=(0,110), yticks=0:10:100)
        plot!(plt, title="Induced Sampling vs Standard Sampling")
        annotate!([1], [100 .* succ[1] ./ trials + 5], text("$(round(100 .* succ[1] ./ trials, digits=1))%"))
        annotate!([2], [100 .* succ[2] ./ trials + 5], text("$(round(100 .* succ[2] ./ trials, digits=1))%"))
    else
        plt = nothing
    end
    return succ, plt
end

function continuous_insertion_trials(d::Union{InducedDistribution,MultivariateInducedDistribution}, trials=100; M0=100, MMax=50M0, tol=1e-6, induced=true)
    allminweights = []
    allconcentrations = []
    allMs = []
    for trial in ProgressBar(1:trials)
        q = least_squares_quad(d, M0, induced=induced)
        minweights = Float64[]
        concentrations = Float64[]
        Ms = Int[]
        for _ in M0:MMax
            add_node!(q)
            push!(minweights, minimum(q.weights))
            push!(concentrations, concentration(q))
            push!(Ms, q.Mtot)
            if is_positive(q, tol=tol)
                break
            end
        end
        push!(allminweights, minweights)
        push!(allconcentrations, concentrations)
        push!(allMs, Ms)
    end
    plt1 = plot()
    plt2 = plot()
    for i in eachindex(allMs)
        plot!(plt1, allMs[i], allminweights[i], alpha=0.5, c=:blue, label=false)
        plot!(plt2, allMs[i], allconcentrations[i], alpha=0.5, c=:red, label=false)
    end
    return plt1, plt2, allminweights, allconcentrations, allMs
end

function continuous_removal_trials(d::Union{InducedDistribution,MultivariateInducedDistribution}, trials=100; M0=100, MMax=50M0, tol=1e-6, induced=true)
    allminweights = []
    allconcentrations = []
    allMs = []
    for trial in ProgressBar(1:trials)
        q = least_squares_quad(d, M0, induced=induced)
        minweights = Float64[]
        concentrations = Float64[]
        Ms = Int[]
        for _ in M0:MMax
            remove_add_node!(q)
            push!(minweights, minimum(q.weights))
            push!(concentrations, concentration(q))
            push!(Ms, q.Mtot)
            if is_positive(q, tol=tol)
                break
            end
        end
        push!(allminweights, minweights)
        push!(allconcentrations, concentrations)
        push!(allMs, Ms)
    end
    plt1 = plot()
    plt2 = plot()
    for i in eachindex(allMs)
        plot!(plt1, allMs[i], allminweights[i], alpha=0.5, c=:blue, label=false)
        plot!(plt2, allMs[i], allconcentrations[i], alpha=0.5, c=:red, label=false)
    end
    return plt1, plt2, allminweights, allconcentrations, allMs
end