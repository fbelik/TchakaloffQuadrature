struct GivensUpDowndateMatrix{T} <: AbstractMatrix{T} 
    V::AbstractMatrix{T}
    Q::AbstractMatrix{T}
    R::AbstractMatrix{T}
end

function GivensUpDowndateMatrix(V::AbstractMatrix)
    M, N = size(V)
    Q,R = qr(V)
    Q = collect(Q)
    R = vcat(R, zeros(M - N, N))
    return GivensUpDowndateMatrix(V,Q,R)
end

Base.size(v::GivensUpDowndateMatrix) = size(v.V)
Base.getindex(v::GivensUpDowndateMatrix, idx::Vararg{Int,2}) = getindex(v.V, idx[1], idx[2])

function reset!(v::GivensUpDowndateMatrix)
    M, N = size(v.V)
    Q,R = qr(v.V)
    v.Q .= collect(Q)
    v.R .= vcat(R, zeros(M - N, N))
    nothing
end

# Copied from CaratheodoryPruning.jl
function givens_qr_row_downdate!(v::GivensUpDowndateMatrix, rowidx)
    Q, R = v.Q, v.R
    n = size(Q,1)
    q = view(Q, rowidx, 1:n)
    r = q[n]
    for i in n:-1:(rowidx+1)
        G, r = givens(q[i-1], r, i-1, i)
        rmul!(Q, G')
        lmul!(G, R)
    end
    for i in (rowidx-1):-1:1
        G, r = givens(r, q[i], rowidx, i)
        rmul!(Q, G')
        lmul!(G, R)
    end
end

# Copied from CaratheodoryPruning.jl
function givens_qr_row_update!(v::GivensUpDowndateMatrix, rowidx, newrow)
    v.V[rowidx,:] .= newrow
    Q, R = v.Q, v.R
    n = size(R,2)
    Q[rowidx, rowidx] = 1.0
    R[rowidx, :] .= newrow
    for i in 1:min(rowidx-1, n)
        G, r = givens(R[i,i], R[rowidx,i], i, rowidx)
        rmul!(Q, G')
        lmul!(G, R)
    end
    for i in rowidx:n
        G, r = givens(R[i,i], R[i+1,i], i, i+1)
        rmul!(Q, G')
        lmul!(G, R)
    end
end