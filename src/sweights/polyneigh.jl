# This file contains functions for creating spatial weights from polygon data

# Auxiliar function to return the points from a vector of polygons, multipolygon, or a mix of the two
function getpointsPoligon(P::Vector)::Tuple{Matrix{Vector{Float64}}, Matrix{Vector{Float64}}}
    n = length(P)
    ppointsx = copy.(fill(Float64[], n, 1))
    ppointsy = copy.(fill(Float64[], n, 1))
        
    for i in 1:n
        
        c = coordinates(P[i])

        ppointsx[i] =  Array{Float64}(undef, length(c[1]))
        ppointsy[i] =  Array{Float64}(undef, length(c[1]))
        
        ctype = typeof(c)

        if ctype == Vector{Vector{Vector{Float64}}}
            # Polygon
            x = map(a -> a[1], c[1])
            y = map(a -> a[2], c[1])

            ppointsx[i] = x
            ppointsy[i] = y
        elseif ctype == Vector{Vector{Vector{Vector{Float64}}}}
            # Multipolygon
            xy = reduce(vcat,
                map(1:length(c)) do i
                    x = map(a -> a[1], c[i][1])
                    y = map(a -> a[2], c[i][1])

                    hcat(x, y)
                end
            )

            x = xy[:,1]
            y = xy[:,2]

            ppointsx[i] = x
            ppointsy[i] = y
        end
    end

    return ppointsx, ppointsy

end

# Auxiliar function to return the number of times the points of two polygons hits
function hits(icoords::Vector{Vector{Float64}}, jcoords::Vector{Vector{Float64}}, critthr::Int64, tol::Float64)::Bool
    nhits = 0

    nicoords = length(icoords[1])
    njcoords = length(jcoords[1])

    for pti in 1:nicoords
        for ptj in 1:njcoords
            @inbounds xd = icoords[1][pti] - jcoords[1][ptj]
            if abs(xd) > tol continue end
            
            @inbounds yd = icoords[2][pti] - jcoords[2][ptj]
            if abs(yd) > tol continue end
            
            dist = hypot(xd, yd)
            if dist <= tol
                nhits = nhits + 1
            end
            
            # Leave if number of hits is larger than the required by the criterion
            if nhits >= critthr
                return true
            end
        end
    end
    
    return false
end

"""
    polyneigh(P, criterion = :Queen)
Build a spatial weights object from a vector of polygons `P`.

# Optional Arguments
- `criterion=:Queen`: neighbour criterion. `:Queen` or `:Rook`.
- `tol=0.0`: tolerance for polygon contiguity.

# Examples
```jldoctest
julia> using SpatialDatasets

julia> guerry = sdataset("Guerry");

julia> W = polyneigh(guerry.geometry)
Spatial Weights object 
Observations: 85 
Transformation: binary
Average number of neighbors: 4.9412
Minimum nunmber of neighbors: 2
Maximum nunmber of neighbors: 8
Median number of neighbors: 5.0
Islands (isloated): 0
Density: 5.8131% 
```
"""
function polyneigh(P::Vector; criterion::Symbol = :Queen, tol::Float64 = 0.0)::SpatialWeights

    n = size(P, 1)

    x, y = getpointsPoligon(P)

    # Bounding box
    xmin, ymin = minimum.(x) .- tol, minimum.(y) .- tol
    xmax, ymax = maximum.(x) .+ tol, maximum.(y) .+ tol

    BBpols = hcat(xmin, ymin, xmax, ymax, [1:n;])

    # Select candidates if bouning box overlaps or touches    
    BBpols = sortslices(BBpols, dims=1)
            
    candidates = copy.(fill(Int[], n, 1))
    for i in 1:n-1
        for j in (i+1):n
            if (BBpols[j,1] <= BBpols[i,3]) && # xmin[j] <= xmax[i]
                (BBpols[j,2] <= BBpols[i,4]) && # ymin[j] <= ymax[i]
                (BBpols[j,3] >= BBpols[i,1]) && # xmax[j] >= xmin[i]
                (BBpols[j,4] >= BBpols[i,2]) # ymax[j] >= ymin[i]
                    push!(candidates[Int(BBpols[i,5])], BBpols[j,5])
            end
        end
    end
    
    # Check neighbours
    neighs = copy.(fill(Int[], n, 1))
    weights = copy.(fill(Float64[], n, 1))
    nneights = zeros(Int,n)

    if criterion == :Queen critthr = 1 end
    if criterion == :Rook  critthr = 2 end

    for i in 1:n
        @inbounds for j in candidates[i]
            p1 = [x[i], y[i]]
            p2 = [x[j], y[j]]
            polhits = hits(p1, p2, critthr, tol)
            
            if polhits  
                push!(neighs[i], j) 
                push!(neighs[j], i) 

                push!(weights[i], 1)
                push!(weights[j], 1)     
                
                nneights[i] += 1
                nneights[j] += 1
              end
        end        
    end

    # Sort neighbours lists
    for i in 1:n
        sort!(neighs[i])
    end

    SpatialWeights(n, neighs, weights, nneights, :binary)
end

