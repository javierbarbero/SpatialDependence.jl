# This file contains functions for creating spatial weights from points data using a distane threshold

"""
    dnearneigh(X, Y; threshold)
Build a spatial weights object from a set of coordinates using a distance `threshold`.

# Examples
```jldoctest
julia> using SpatialDatasets

julia> bostonhsg = sdataset("Bostonhsg");

julia> W = dnearneigh(bostonhsg.x, bostonhsg.y, threshold = 4)
Spatial Weights object 
Observations: 506 
Transformation: binary
Average number of neighbors: 72.2253
Minimum nunmber of neighbors: 1
Maximum nunmber of neighbors: 178
Median number of neighbors: 56.0
Islands (isloated): 0
Density: 14.2738% 
```
"""
function dnearneigh(X::Vector{Float64}, Y::Vector{Float64}; threshold)::SpatialWeights

    n, ny = length(X), length(Y)

    n == ny || throw(DimensionMismatch("dimensions must match: X has length ($(n)), Y has length ($ny)"))
    
    # Get points within range using a KDTree
    cpoints = vcat(X', Y')

    kdtree = KDTree(cpoints)    
    idxs  = inrange(kdtree, cpoints, threshold, true)
        
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
        
    SpatialWeights(n, neighs, weights, nneights, :binary)

end

"""
    dnearneigh(P; threshold)
Build a spatial weights object from a vector of points using a distance `threshold`.

# Examples
```jldoctest
julia> using SpatialDatasets

julia> bostonhsg = sdataset("Bostonhsg");

julia> W = dnearneigh(bostonhsg.geometry, threshold = 4)
Spatial Weights object 
Observations: 506 
Transformation: binary
Average number of neighbors: 72.2253
Minimum nunmber of neighbors: 1
Maximum nunmber of neighbors: 178
Median number of neighbors: 56.0
Islands (isloated): 0
Density: 14.2738% 
```
"""
function dnearneigh(P::Vector; threshold)::SpatialWeights

    cpoints = reduce(hcat, map(a -> coordinates(a), P))

    dnearneigh(cpoints[1,:], cpoints[2,:]; threshold)
end
