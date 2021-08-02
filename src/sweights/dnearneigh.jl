# This file contains functions for creating spatial weights from points data using a distane threshold

"""
    dnearneigh(X, Y; threshold)
Build a spatial weights object from a set of coordinates using a distance `threshold`.
"""
function dnearneigh(X::Vector{Float64}, Y::Vector{Float64}; threshold::Real)::SpatialWeights
    n, ny = length(X), length(Y)

    n == ny || throw(DimensionMismatch("dimensions must match: X has length ($(n)), Y has length ($ny)"))
    
    # Get points within range using a KDTree
    cpoints = vcat(X', Y')

    kdtree = KDTree(cpoints)    
    idxs  = inrange(kdtree, cpoints, threshold, true)
        
    # Check neighbours
    neighs = copy.(fill(Int[], n))
    weights = copy.(fill(Float64[], n))
    nneighs = zeros(Int,n)
        
    for i in 1:n   
        neighsi = idxs[i]    
        neighs[i] = sort!(neighsi[neighsi .!= i])
        nneighs[i] = length(neighs[i])
        weights[i] = ones(nneighs[i]) ./ nneighs[i]
    end
        
    SpatialWeights(n, neighs, weights, nneighs, :row)

end

"""
    dnearneigh(P; threshold)
Build a spatial weights object from a vector of points using a distance `threshold`.
"""
function dnearneigh(P::Vector{T} where T<:Union{Missing,AbstractPoint}; threshold::Real)::SpatialWeights
    cpoints = reduce(hcat, map(a -> coordinates(a), P))

    dnearneigh(cpoints[1,:], cpoints[2,:]; threshold)
end

"""
    dnearneigh(A; threshold)
Build a spatial weights object from a table A that constains a points geometry column using a distance `threshold`.
"""
function dnearneigh(A::Any; threshold::Real)::SpatialWeights
    istable(A) || throw(ArgumentError("Unknown geometry or not points geometry"))

    (:geometry in propertynames(A)) || throw(ArgumentError("table does not have :geometry information"))

    dnearneigh(A.geometry; threshold)
end
