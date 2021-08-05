# This file contains the structure and common functions for SpatialWeights
mutable struct SpatialWeights
    n::Int64
    neighs::Vector{Vector{Int64}}
    weights::Vector{Vector{Float64}}
    nneighs::Vector{Int64}
    transform::Symbol
end

function SpatialWeights(W::AbstractMatrix; standardize = true)

    n, nj = size(W)

    n == nj || throw(ArgumentError("matrix is not square"))
    
    neighs = copy.(fill(Int[], n))
    weights = copy.(fill(Float64[], n))
    nneighs = zeros(Int, n)

    for i in 1:n
        neighs[i] = (1:n)[W[i,:] .!= 0]
        nneighs[i] = length(neighs[i])
        if standardize
            weights[i] = ones(nneighs[i]) ./ nneighs[i]
        else
            weights[i] = W[i,neighs[i]]
        end        
    end

    SpatialWeights(n, neighs, weights, nneighs, standardize ? :row : :original)

end

"""
    cardinalities(W::SpatialWeights)
Return a vector with the number of neighbors of each observation.
"""
cardinalities(W::SpatialWeights)::Vector{Int64} = W.nneighs;

"""
    neighbors(W::SpatialWeights, i::Int)
Return a vector with the neighbors of ``i``.
"""
neighbors(W::SpatialWeights, i::Int)::Vector{Int64} = W.neighs[i];

"""
    weights(W::SpatialWeights, i::Int)
Return a vector with the weights of the neighbors of ``i``.
"""
weights(W::SpatialWeights, i::Int)::Vector{Float64} = W.weights[i];

"""
    wtransform!(W::SpatialWeights, style::Symbol)
In-place transformation of the weights using the specified ``style``.

# Weights transformation
`style` can be one of the following:
- `:binary`: 1 if neighbor, 0 if not.
- `:row`: row-standardized. Each row sum equals one.
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
"""
wtransformation(W::SpatialWeights)::Symbol = W.transform ;

"""
    nislands (W::SpatialWeights)
Return the number of islands in the spatial weights object.
"""
nislands(W::SpatialWeights)::Int64 = sum(W.nneighs .== 0) ;

"""
    islands (W::SpatialWeights)
Return a vector with the islands in the spatial weights object.
"""
function islands(W::SpatialWeights)::Vector{Int64}
    islandsvec = Int[]
    for i = 1:W.n
        if W.nneighs[i] == 0
            push!(islandsvec, i)
        end
    end

    return islandsvec
end

# Descriptive statistics
nobs(W::SpatialWeights)::Int64 = W.n;

mean(W::SpatialWeights)::Float64 = mean(W.nneighs);

Base.minimum(W::SpatialWeights)::Int64 = minimum(W.nneighs);

Base.maximum(W::SpatialWeights)::Int64 = maximum(W.nneighs);

median(W::SpatialWeights)::Float64 = median(W.nneighs);

# Transformation functions
function Base.Matrix(W::SpatialWeights)::Matrix{Float64}
    n = W.n

    mat = zeros(n, n)
    for i in 1:n        
        mat[i, W.neighs[i]] = W.weights[i]
    end

    return mat

end

function sparse(W::SpatialWeights)::SparseMatrixCSC{Float64, Int64}
    n = W.n
    tn = sum(W.nneighs)

    # Build I, J and V vectors to build the spatial matrix
    I = zeros(Int, tn)
    J = zeros(Int, tn)
    V = zeros(tn)

    nlast = 0
    for i in 1:n
        sub = (nlast + 1):(nlast + W.nneighs[i])
        I[sub] = fill(i, W.nneighs[i])
        J[sub] = W.neighs[i]
        V[sub] = W.weights[i]

        nlast = nlast + W.nneighs[i]
    end

    return sparse(I, J, V)

end

function Base.show(io::IO, W::SpatialWeights)

    print(io, "Spatial Weights \n")
    print(io, "Observations: ", nobs(W), " \n")    
    print(io, "Transformation: ", wtransformation(W), "\n")
    print(io, "Minimum nunmber of neighbors: ", minimum(W), "\n")
    print(io, "Maximum nunmber of neighbors: ", maximum(W), "\n")
    print(io, "Average number of neighbors: ", round(mean(W), digits = 4), "\n")
    print(io, "Median number of neighbors: ", median(W), "\n")
    print(io, "Islands (isloated): ", nislands(W), "\n")
    print(io, "Density: ", round(100 * sum(W.nneighs) / (W.n  * W.n), digits = 4), "% \n")

end
