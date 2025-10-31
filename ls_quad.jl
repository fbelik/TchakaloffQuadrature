mutable struct LS_Quad{T}
    Mtot::Int
    Mcur::Int
    const d::Union{InducedDistribution,MultivariateInducedDistribution}
    V::Matrix{T}
    const square_pairs::Vector{Tuple{Int,Int}}
    Vsq::AbstractMatrix{T}
    const wsq::Vector{T}
    const τ²::Vector{T}
    const CS_rescale::Vector{T}
    const Gm::Matrix{T}
    Gminv::Union{Matrix{T},Nothing}
    const η::Vector{T}
    const nodes::Vector
    const weights::Vector{T}
    const induced::Bool
end

# Assuming G=I
function least_squares_quad(d::InducedDistribution, M=nothing; induced=true, poly_square_set=true)
    N = d.N+1
    square_pairs = begin
        if poly_square_set
            index_set_square_pairs([(i,) for i in 0:d.N])
        else
            [(i,j) for i in 1:N for j in i:N]
        end
    end
    if isnothing(M)
        M = length(square_pairs) + 1
    end
    if induced
        nodes = rand(d, M)
    else
        d0 = _dists[d.poly]
        nodes = rand(d0, M)
    end
    V = evaluate(0:d.N, nodes, d.PolyObj) ./ sqrt.(d.normssq)'
    τ² = sum(V .^ 2, dims=2)[:,1] ./ N
    if !induced
        τ² .= 1.0
    end
    CS_rescale = ones(M)
    Vsq = Matrix{Float64}(undef, M, length(square_pairs))
    for (idx,(i,j)) in enumerate(square_pairs)
        for m in 1:M
            Vsq[m,idx] = V[m,i] * V[m,j] / τ²[m]
        end
    end
    wsq = ones(M) ./ M
    Gm = (V' * Diagonal(1 ./ τ²) * V) ./ M
    η = d.V * d.w_Q
    w = zeros(M)
    Gminv = nothing
    try
        w .= transpose((transpose(η) / Gm) * V' * Diagonal(CS_rescale ./ τ²)) ./ M
        Gminv = inv(Gm)
    catch e
        if !isa(e, SingularException)
            throw(e)
        end
        w .= transpose((transpose(η) / pinv(Gm)) * V' * Diagonal(CS_rescale ./ τ²)) ./ M
    end
    return LS_Quad(M, M, d, V, square_pairs, Vsq, wsq, τ², CS_rescale, Gm, Gminv, η, nodes, w, induced)
end

function least_squares_quad(d::MultivariateInducedDistribution, M=nothing; induced=true, poly_square_set=true)
    N = length(d.index_set)
    square_pairs = begin
        if poly_square_set
            index_set_square_pairs(d.index_set)
        else
            [(i,j) for i in 1:N for j in i:N]
        end
    end
    if isnothing(M)
        M = length(square_pairs) + 1
    end
    D = length(d.polys)
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
    CS_rescale = ones(M)
    Vsq = Matrix{Float64}(undef, M, length(square_pairs))
    for (idx,(i,j)) in enumerate(square_pairs)
        for m in 1:M
            Vsq[m,idx] = V[m,i] * V[m,j] / τ²[m]
        end
    end
    wsq = ones(M) ./ M
    η = zeros(N); η[1] = 1
    nodes = [collect(x) for x in eachcol(nodes)]
    Gm = (V' * Diagonal(1 ./ τ²) * V) ./ (M)
    G = I
    w = transpose((transpose(η) / Gm) * V' * Diagonal(CS_rescale ./ τ²)) ./ M
    Gminv = nothing
    try
        Gminv = inv(Gm)
    catch e
        if !isa(e, SingularException)
            throw(e)
        end
    end
    return LS_Quad(M, M, d, V, square_pairs, Vsq, wsq, τ², CS_rescale, Gm, Gminv, η, nodes, w, induced)
end

function add_node!(q::LS_Quad, x=nothing; nodeidx=nothing, gminvtol=1e-2)
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
    if q.Mcur < q.Mtot
        # Rescale CS_rescale
        q.CS_rescale .*= (q.Mcur + 1) / q.Mcur
        q.CS_rescale .*= q.Mtot / (q.Mtot + 1)
    end
    if isnothing(nodeidx)
        # Add new row to V
        q.V = vcat(q.V, v')
        # Update τ²
        push!(q.τ², τ²new)
        # Update CS_rescale
        push!(q.CS_rescale, (q.Mcur+1)/(q.Mtot+1))
        # Add new row to Vsq
        q.Vsq = vcat(q.Vsq, Vector{Float64}(undef, length(q.square_pairs))')
        for (idx,(i,j)) in enumerate(q.square_pairs)
            q.Vsq[end,idx] = v[i] * v[j] / τ²new
            idx += 1
        end
        # Add new element to wsq
        q.wsq .*= q.Mtot / (q.Mtot+1)
        push!(q.wsq, 1 / (q.Mtot+1))
        # Add node
        push!(q.nodes, x)
        # Insert new weight
        push!(q.weights, 0)
    else
        # Insert new row to V
        q.V[nodeidx,:] .= v
        # Update τ²
        q.τ²[nodeidx] = τ²new
        # Update CS_rescale
        q.CS_rescale[nodeidx] = (q.Mcur+1)/(q.Mtot+1)
        # Insert new row to Vsq
        if !isa(q.Vsq, GivensUpDowndateMatrix)
            q.Vsq = GivensUpDowndateMatrix(q.Vsq)
        end
        newrow = Vector{Float64}(undef, size(q.Vsq, 2))
        for (idx,(i,j)) in enumerate(q.square_pairs)
            newrow[idx] = v[i] * v[j] / τ²new
        end
        givens_qr_row_update!(q.Vsq, nodeidx, newrow)
        # Update wsq
        q.wsq .*= q.Mtot / (q.Mtot+1)
        q.wsq[nodeidx] = 1 / (q.Mtot+1)
        # Insert node
        q.nodes[nodeidx] = x
    end
    # Update Gm
    q.Gm .= q.Mtot / (q.Mtot+1) .* q.Gm .+ (v * v') ./ (τ²new * (q.Mtot+1))
    # Update Gm⁻¹ by Sherman-Morrison
    if !isnothing(q.Gminv)
        prod = q.Gminv * v
        q.Gminv .= (q.Mtot+1)/q.Mtot .* (q.Gminv .- (prod * prod') ./ (q.Mtot*τ²new + dot(v,q.Gminv,v)))
    else
        try
            q.Gminv = inv(q.Gm)
        catch e
            if !isa(e, SingularException)
                throw(e)
            end
            q.Gminv = nothing
        end
    end
    # Update M
    q.Mcur += 1
    q.Mtot += 1
    # Update weights
    q.weights .= transpose((transpose(q.η) * q.Gminv) * q.V' * Diagonal(q.CS_rescale ./ q.τ²)) ./ q.Mcur
    # Check Gminv
    if !isnothing(q.Gminv) && norm(q.Gm * q.Gminv - I) > gminvtol
        try
            q.Gminv .= inv(q.Gm)
        catch e
            if !isa(e, SingularException)
                throw(e)
            end
            q.Gminv = nothing
        end
    end
    return
end

function remove_node!(q::LS_Quad; gminvtol=1e-2, fullqr=false, qrresetpct=1e-3)
    if size(q.Vsq, 1) > size(q.Vsq, 2) && !isnothing(q.Gminv)
        # Remove node and resize
        # Select node to remove from q
        if !isa(q.Vsq, GivensUpDowndateMatrix)
            q.Vsq = GivensUpDowndateMatrix(q.Vsq)
        elseif fullqr # Recompute QR
            q.Vsq = GivensUpDowndateMatrix(q.Vsq.V)
        end
        if rand() < qrresetpct
            reset!(q.Vsq)
        end
        n = q.Vsq.Q[:,end]
        _, remidx = findmin(abs.(q.wsq ./ n))
        α = q.wsq[remidx] / n[remidx]
        q.wsq .-= α .* n
        # Downdate q.Vsq
        givens_qr_row_downdate!(q.Vsq, remidx)
        # Remove node remidx from q
        new_indices = [1:remidx-1 ; remidx+1:size(q.V,1)]
        # Update weights, and CS_rescale
        q.CS_rescale .*= (q.wsq ./ (q.wsq .+ α .* n))
        q.Mcur -= 1
        q.CS_rescale .*= q.Mcur / (q.Mcur + 1)
        q.weights[new_indices] .= transpose((transpose(q.η) * q.Gminv) * view(q.V,new_indices,:)' * Diagonal(view(q.CS_rescale, new_indices) ./ view(q.τ²,new_indices))) ./ q.Mcur
        # Resize arrays
        q.V = q.V[new_indices,:]
        q.Vsq = q.Vsq[new_indices,:]
        deleteat!(q.wsq, remidx)
        deleteat!(q.τ², remidx)
        deleteat!(q.CS_rescale, remidx)
        deleteat!(q.nodes, remidx)
        deleteat!(q.weights, remidx)
        return
    end
end

function remove_add_node!(q::LS_Quad, x=nothing; gminvtol=1e-2, fullqr=false, qrresetpct=1e-3)
    if size(q.Vsq, 1) > size(q.Vsq, 2) && !isnothing(q.Gminv)
        # Remove and add
        # Select node to remove from q
        if !isa(q.Vsq, GivensUpDowndateMatrix)
            q.Vsq = GivensUpDowndateMatrix(q.Vsq)
        elseif fullqr # Recompute QR
            q.Vsq = GivensUpDowndateMatrix(q.Vsq.V)
        end
        if rand() < qrresetpct
            reset!(q.Vsq)
        end
        n = q.Vsq.Q[:,end]
        _, remidx = findmin(abs.(q.wsq ./ n))
        α = q.wsq[remidx] / n[remidx]
        q.wsq .-= α .* n
        # Downdate q.Vsq
        givens_qr_row_downdate!(q.Vsq, remidx)
        # Remove node remidx from q
        new_indices = [1:remidx-1 ; remidx+1:size(q.V,1)]
        # Update weights, and CS_rescale
        q.CS_rescale .*= (q.wsq ./ (q.wsq .+ α .* n))
        q.Mcur -= 1
        q.CS_rescale .*= q.Mcur / (q.Mcur + 1)
        q.weights[new_indices] .= transpose((transpose(q.η) * q.Gminv) * view(q.V,new_indices,:)' * Diagonal(view(q.CS_rescale, new_indices) ./ view(q.τ²,new_indices))) ./ q.Mcur
        # Add new node at remidx
        add_node!(q, x, nodeidx=remidx, gminvtol=gminvtol)
    else
        # Simply add a node
        add_node!(q, x, gminvtol=gminvtol)
    end
end

function moment_error(q::LS_Quad, norm=norm)
    return norm(transpose(q.V) * q.weights .- q.η)
end

function ls_nodes_weights_pruned(q::LS_Quad; kwargs...)
    V = q.V
    w = q.weights
    w_pruned,inds = caratheodory_pruning(V, w; kwargs...) 
    return q.nodes[inds], w_pruned[inds], inds
end

function ls_weights(q::LS_Quad, CS_rescaled=true)
    if CS_rescaled
        return q.weights
    else
        Gm = (q.V' * Diagonal(1 ./ q.τ²) * q.V) ./ M
        return transpose((transpose(q.η) / Gm) * q.V' * Diagonal(1 ./ q.τ²)) ./ q.Mcur
    end
end

function monte_carlo_weights(q::LS_Quad, CS_rescaled=true)
    if CS_rescaled
        return q.wsq ./ q.τ²
    else
        return 1 / q.Mcur ./ q.τ²
    end
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

function is_positive(q::LS_Quad, CS_rescaled=true; tol=0.0)
    mval = -1 * abs(tol)
    for w in ls_weights(q, CS_rescaled)
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

function continuous_insertion_trials(d::Union{InducedDistribution,MultivariateInducedDistribution}, trials=100; MMax=5000, tol=1e-6, induced=true)
    allminweights = []
    allconcentrations = []
    allMs = []
    for trial in ProgressBar(1:trials)
        q = least_squares_quad(d, induced=induced)
        minweights = Float64[]
        concentrations = Float64[]
        Ms = Int[]
        while q.Mtot < MMax
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

function continuous_removal_trials(d::Union{InducedDistribution,MultivariateInducedDistribution}, trials=100; MMax=5000, tol=1e-6, induced=true)
    allminweights = []
    allconcentrations = []
    allMs = []
    for trial in ProgressBar(1:trials)
        q = least_squares_quad(d, induced=induced)
        minweights = Float64[]
        concentrations = Float64[]
        Ms = Int[]
        while q.Mtot < MMax
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