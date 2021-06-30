# This file contains functions for calculating spatial lags

"""
    slag(W, x)
Calculate the spatial lag of `x` using the spatial weights `W`.

# Examples
```jldoctest
julia> using SpatialDatasets

julia> guerry = sdataset("Guerry");

julia> W = polyneigh(guerry.geometry);

julia> wtransform!(W, :row)

julia> sx = slag(W, guerry.Litercy);
```
"""
function slag(W::SpatialWeights, x::Vector{Float64})::Vector{Float64}
    n = W.n
    nx = length(x)

    n == nx || throw(DimensionMismatch("dimensions must match: W has ($(n)), x has ($nx)"))
        
    sx = zeros(n)    
    for i in 1:n
        @inbounds ni = length(W.neighs[i])
        for j in 1:ni
            @inbounds sx[i] = sx[i] + W.weights[i][j] * x[W.neighs[i][j]]
        end
    end
    
    return sx

end