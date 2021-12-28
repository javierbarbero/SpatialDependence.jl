# This file contains functions for calculating spatial lags

"""
    slag(W, x)
Calculate the spatial lag of `x` using the spatial weights `W`.
"""
function slag(W::SpatialWeights, x::AbstractVector{T}) where {T <:Real}
    n = W.n
    nx = length(x)

    n == nx || throw(DimensionMismatch("dimensions must match: W has ($(n)), x has ($nx)"))
        
    sx = zeros(n)    
    @inbounds @simd for i in 1:n
        ni = length(W.neighs[i])
        for j in 1:ni
            sx[i] = sx[i] + W.weights[i][j] * x[W.neighs[i][j]]
        end
    end
    
    return sx

end

Base.:*(W::SpatialWeights, x::AbstractVector{T}) where {T <:Real} = slag(W, x);
