function in_td_set(idx, N=5)
    return sum(idx) <= N
end

function in_tp_set(idx, N=5)
    return maximum(idx) <= N
end

function in_hc_set(idx, N=5)
    return prod(idx .+ 1) <= (N + 1)
end

in_set_dict = Dict(
    :hc => in_hc_set,
    :hyperbolic_cross => in_hc_set,
    :td => in_td_set,
    :total_degree => in_td_set,
    :tp => in_tp_set,
    :total_product => in_tp_set
)

function multi_index_set(D=2,N=5,set=:td)
    in_set = in_set_dict[set]
    idx = zeros(Int, D)
    index_set = [Tuple(idx)]
    idx[1] += 1
    while true
        # Add to index set
        push!(index_set, Tuple(idx))
        # Update idx
        for i in 1:D
            idx[i] += 1
            if in_set(idx, N)
                break
            end
            idx[i] = 0
        end
        if all(iszero, idx)
            break
        end
    end
    return index_set
end

function visualize_multi_index_set(set; kwargs...)
    xs = [idx[1] for idx in set]
    ys = [idx[2] for idx in set]
    if length(set[1]) > 2
        zs = [idx[3] for idx in set]
        return scatter(xs, ys, zs, label=false, title="Multi Index Set"; kwargs...)
    end
    return scatter(xs, ys, label=false, title="Multi Index Set"; kwargs...)
end