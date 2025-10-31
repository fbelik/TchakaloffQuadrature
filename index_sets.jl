function TD(idx, N=5)
    return sum(idx) <= N
end

function TP(idx, N=5)
    return maximum(idx) <= N
end

function HC(idx, N=5)
    return prod(idx .+ 1) <= (N + 1)
end

function multi_index_set(D=2,N=5,in_set=TD)
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

function index_set_square_pairs(index_set::Vector{<:Tuple})
    index_set_squared = [(i,j,i1 .+ i2) for (i,i1) in enumerate(index_set) for (j,i2) in enumerate(index_set)]
    unique!(x -> x[3], index_set_squared) # Remove duplicates
    square_pairs = [x[1:2] for x in index_set_squared]
    return square_pairs
end

function square(index_set::Vector{<:Tuple})
    index_set_squared = [i1 .+ i2 for i1 in index_set for i2 in index_set]
    return unique(index_set_squared)
end

function multi_index_set_squared(D=2,N=5,in_set=TD;exact=true)
    if in_set in (TP, TD)
        return multi_index_set(D, 2N, in_set)
    elseif in_set == HC
        if exact
            return square(multi_index_set(D, N, in_set))
        else
            return multi_index_set(D, (N+1)^2-1, in_set)
        end
    else
        # Unknown, do exact
        return square(multi_index_set(D, N, in_set))
    end
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

function visualize_multi_index_set!(set; kwargs...)
    xs = [idx[1] for idx in set]
    ys = [idx[2] for idx in set]
    if length(set[1]) > 2
        zs = [idx[3] for idx in set]
        return scatter!(xs, ys, zs, label=false, title="Multi Index Set"; kwargs...)
    end
    return scatter!(xs, ys, label=false, title="Multi Index Set"; kwargs...)
end