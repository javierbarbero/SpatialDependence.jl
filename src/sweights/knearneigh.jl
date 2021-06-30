# This file contains functions for creating spatial weights from points data using a k nearest neighbors

"""
    knearneigh(X, Y; k)
Build a spatial weights object from a set of coordinates using `k` nearest neighbors.

# Examples
```jldoctest
julia> using SpatialDatasets

julia> bostonhsg = sdataset("Bostonhsg");

julia> W = knearneigh(bostonhsg.x, bostonhsg.y, k = 10)
Spatial Weights object 
Observations: 506 
Transformation: Binary
Average number of neighbors: 10.0
Minimum nunmber of neighbors: 10
Maximum nunmber of neighbors: 10
Median number of neighbors: 10.0
Islands (isloated): 0
Density: 1.9763% 
```
"""
function knearneigh(X::Vector{Float64}, Y::Vector{Float64}; k)::SpatialWeights

    n, ny = length(X), length(Y)

    n == ny || throw(DimensionMismatch("dimensions must match: X has length ($(n)), Y has length ($ny)"))
    
    # Get nearest neighbors using a KDTree
    cpoints = vcat(X', Y')

    kdtree = KDTree(cpoints)    
    idxs,  = knn(kdtree, cpoints, k + 1, true)
        
    # Check neighbours
    neighs = copy.(fill(Int[], n, 1))
    weights = copy.(fill(Float64[], n, 1))
    nneights = zeros(Int,n)
        
    for i in 1:n   
        neighsi = idxs[i]    
        neighs[i] = sort!(neighsi[neighsi .!= i])
        nneights[i] = length(neighs[i])
        weights[i] = ones(nneights[i])        
    end
        
    SpatialWeights(n, neighs, weights, nneights, :Binary)

end

"""
    knearneigh(P; k)
Build a spatial weights object from a vector of points using `k` nearest neighbors.

# Examples
```jldoctest
julia> using SpatialDatasets

julia> bostonhsg = sdataset("Bostonhsg");

julia> W = knearneigh(bostonhsg.geometry, k = 10)
Spatial Weights object 
Observations: 506 
Transformation: Binary
Average number of neighbors: 10.0
Minimum nunmber of neighbors: 10
Maximum nunmber of neighbors: 10
Median number of neighbors: 10.0
Islands (isloated): 0
Density: 1.9763% 
```
"""
function knearneigh(P::Vector; k)::SpatialWeights

    cpoints = reduce(hcat, map(a -> coordinates(a), P))

    knearneigh(cpoints[1,:], cpoints[2,:]; k)
end
