# This file contains functions for creating spatial weights from points data using a k nearest neighbors

"""
    knearneigh(X, Y; k)
Build a spatial weights object from a set of coordinates using `k` nearest neighbors.
"""
function knearneigh(X::Vector, Y::Vector; k::Int)::SpatialWeights
    n, ny = length(X), length(Y)

    n == ny || throw(DimensionMismatch("dimensions must match: X has length ($(n)), Y has length ($ny)"))
    
    # Fail if missing values
    !any(ismissing.(X)) || throw(DimensionMismatch("missing values not allowed in X coordinates"))
    !any(ismissing.(Y)) || throw(DimensionMismatch("missing values not allowed in Y coordinates"))

    # Get nearest neighbors using a KDTree
    cpoints = vcat(collect(skipmissing(X))', collect(skipmissing(Y))')

    kdtree = KDTree(cpoints)    
    idxs,  = knn(kdtree, cpoints, k + 1, true)
        
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
    knearneigh(P; k)
Build a spatial weights object from a vector of points using ``k`` nearest neighbors.
"""
function knearneigh(P::Vector; k::Int)::SpatialWeights
    all(GI.isgeometry.(P)) || throw(ArgumentError("Unknown geometry"))
    all(isa.(GI.geomtrait.(P), GI.PointTrait)) || throw(ArgumentError("Geometry must be PointTrait"))

    cpoints = reduce(hcat, map(a -> GI.coordinates(a), P))

    knearneigh(cpoints[1,:], cpoints[2,:]; k)
end

"""
    knearneigh(A; k)
Build a spatial weights object from a table A that constains a points geometry column using ``k`` nearest neighbors..
"""
function knearneigh(A::Any; k::Int)::SpatialWeights
    istable(A) || throw(ArgumentError("Argument must be a table with geometry or a vector of points"))

    geomcol = _geomFromTable(A)

    knearneigh(geomcol; k)
end
