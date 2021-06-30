# This file contains the structure and common functions for SpatialWeights
mutable struct SpatialWeights
    n::Int64
    neighs::Matrix{Vector{Int64}}
    weights::Matrix{Vector{Float64}}
    nneights::Vector{Int64}
    transform::Symbol
end

function SpatialWeights(W::AbstractMatrix)

    n, nj = size(W)

    n == nj || throw(ArgumentError("matrix is not square"))
    
    neighs = copy.(fill(Int[], n, 1))
    weights = copy.(fill(Float64[], n, 1))
    nneights = zeros(Int, n)

    for i in 1:n
        neighs[i] = (1:n)[W[i,:] .!= 0]
        weights[i] = W[i,neighs[i]]
        nneights[i] = length(neighs[i])
    end

    SpatialWeights(n, neighs, weights, nneights, :original)

end

"""
    neighbors(W::SpatialWeights, i::Int)
Return a vector with the neighbors of `i`.
# Examples
```jldoctest
julia> W = SpatialWeights([0 1 0; 1 0 1; 0 1 0]);

julia> neighbors(W, 2)
2-element Vector{Int64}:
 1
 3
```
"""
neighbors(W::SpatialWeights, i::Int)::Vector{Int64} = W.neighs[i];

"""
    weights(W::SpatialWeights, i::Int)
Return a vector with the weights of the neighbors of `i`.
# Examples
```jldoctest
julia> W = SpatialWeights([0 1 0; 1 0 1; 0 1 0]);

julia> wtransform!(W, :row)

julia> weights(W, 2)
2-element Vector{Float64}:
 0.5
 0.5
```
"""
weights(W::SpatialWeights, i::Int)::Vector{Float64} = W.weights[i];

"""
    wtransform!(W::SpatialWeights, style::Symbol)
In-place transformation of the weights using the specified `style`.

# Weights transformation
`style` can be one of the following:
- `:binary`: 1 if neighbor, 0 if not.
- `:row`: row-standardized. Each row sum equals one.

# Examples
```jldoctest
julia> W = SpatialWeights([0 1 0; 1 0 1; 0 1 0]);

julia> wtransform!(W, :row)
```
"""
function wtransform!(W::SpatialWeights, style::Symbol)

    n = W.n
    
    (style == :binary || style == :row) ||  throw(ArgumentError("unkown transformation $(style)"))

    if style == :binary
        # Binary transformation
        for i in 1:n
            W.weights[i] = ones(length(W.neighs[i]))
        end
    
    else style == :row
        # Row transform
        for i in 1:n
            W.weights[i] = ones(length(W.neighs[i])) ./ length(W.neighs[i])
        end
 
    end

    W.transform = style
    return nothing

end

"""
    wtransform(W::SpatialWeights, style::Symbol)
Returns a transformed copy of the weights matrix.
"""
function wtransform(W::SpatialWeights, style::Symbol)::SpatialWeights
    Wnew = deepcopy(W)
    wtransform!(Wnew, style)
    return Wnew
end

"""
    wtransformation (W::SpatialWeights)
Return the current transformation of the spatial weights.

# Examples
```jldoctest
julia> W = SpatialWeights([0 1 0; 1 0 1; 0 1 0]);

julia> wtransform!(W, :row)

julia> wtransformation(W)
:row
```
"""
wtransformation(W::SpatialWeights)::Symbol = W.transform ;

"""
    nislands (W::SpatialWeights)
Return the number of islands in the spatial weights object.

# Examples
```jldoctest
julia> W = SpatialWeights([0 1 0; 1 0 1; 0 0 0]);

julia> nislands(W)
1
```
"""
nislands(W::SpatialWeights)::Int64 = sum(W.nneights .== 0) ;

"""
    islands (W::SpatialWeights)
Return a vector with the islands in the spatial weights object.

# Examples
```jldoctest
julia> W = SpatialWeights([0 1 0; 1 0 1; 0 0 0]);

julia> islands(W)
1-element Vector{Int64}:
 3
```
"""
function islands(W::SpatialWeights)::Vector{Int64}
    islandsvec = Int[]
    for i = 1:W.n
        if W.nneights[i] == 0
            push!(islandsvec, i)
        end
    end

    return islandsvec
end

# Descriptive statistics
nobs(W::SpatialWeights)::Int64 = W.n;

mean(W::SpatialWeights)::Float64 = mean(W.nneights);

Base.minimum(W::SpatialWeights)::Int64 = minimum(W.nneights);

Base.maximum(W::SpatialWeights)::Int64 = maximum(W.nneights);

median(W::SpatialWeights)::Float64 = median(W.nneights);

# Transformation functions
function Matrix(W::SpatialWeights)::Matrix{Float64}
    n = W.n

    mat = zeros(n, n)
    for i in 1:n        
        mat[i, W.neighs[i]] = W.weights[i]
    end

    return mat

end

function sparse(W::SpatialWeights)::SparseMatrixCSC{Float64, Int64}
    n = W.n
    tn = sum(W.nneights)

    # Build I, J and V vectors to build the spatial matrix
    I = zeros(Int, tn)
    J = zeros(Int, tn)
    V = zeros(tn)

    nlast = 0
    for i in 1:n
        sub = (nlast + 1):(nlast + W.nneights[i])
        I[sub] = fill(i, W.nneights[i])
        J[sub] = W.neighs[i]
        V[sub] = W.weights[i]

        nlast = nlast + W.nneights[i]
    end

    return sparse(I, J, V)

end

function Base.show(io::IO, W::SpatialWeights)

    print(io, "Spatial Weights object \n")
    print(io, "Observations: ", nobs(W), " \n")    
    print(io, "Transformation: ", wtransformation(W), "\n")
    print(io, "Average number of neighbors: ", round(mean(W), digits = 4), "\n")
    print(io, "Minimum nunmber of neighbors: ", minimum(W), "\n")
    print(io, "Maximum nunmber of neighbors: ", maximum(W), "\n")
    print(io, "Median number of neighbors: ", median(W), "\n")
    print(io, "Islands (isloated): ", nislands(W), "\n")
    print(io, "Density: ", round(100 * sum(W.nneights) / (W.n  * W.n), digits = 4), "% \n")

end
