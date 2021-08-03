# This file contains functions for creating spatial weights from polygon data

# Auxiliar function to return the points from a vector of polygons, multipolygon, or a mix of the two
function _getpointsPoligon(P::T where T <:AbstractPolygon)::Tuple{Vector{Float64}, Vector{Float64}}
    c = coordinates(P)

    ppointsx =  Array{Float64}(undef, length(c[1]))
    ppointsy =  Array{Float64}(undef, length(c[1]))

    ppointsx = map(a -> a[1], c[1])
    ppointsy = map(a -> a[2], c[1])

    return ppointsx, ppointsy
end

function _getpointsPoligon(P::T where T <:AbstractMultiPolygon)::Tuple{Vector{Float64}, Vector{Float64}}
    c = coordinates(P)

    ppointsx =  Array{Float64}(undef, length(c[1]))
    ppointsy =  Array{Float64}(undef, length(c[1]))

    xy = reduce(vcat,
                map(1:length(c)) do i
                    x = map(a -> a[1], c[i][1])
                    y = map(a -> a[2], c[i][1])

                    hcat(x, y)
                end
            )

    ppointsx = xy[:,1]
    ppointsy = xy[:,2]

    return ppointsx, ppointsy
end

_getpointsPoligon(::Missing) = throw(ArgumentError("Missing geometry"));

# Auxiliar function to return the number of times the points of two polygons hits
function _hits(icoords::Vector{Vector{Float64}}, jcoords::Vector{Vector{Float64}}, critthr::Int64, tol::Float64)::Bool
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
"""
function polyneigh(P::Vector{T} where T <:Union{Missing,AbstractPolygon,AbstractMultiPolygon}; criterion::Symbol = :Queen, tol::Float64 = 0.0)::SpatialWeights
    n = length(P)

    xy = _getpointsPoligon.(P)

    x = map(a -> a[1], xy)
    y = map(a -> a[2], xy)

    # Bounding box
    xmin, ymin = minimum.(x) .- tol, minimum.(y) .- tol
    xmax, ymax = maximum.(x) .+ tol, maximum.(y) .+ tol

    BBpols = hcat(xmin, ymin, xmax, ymax, [1:n;])

    # Select candidates if bouning box overlaps or touches    
    BBpols = sortslices(BBpols, dims=1)
            
    candidates = copy.(fill(Int[], n))
    Threads.@threads for i in 1:n-1
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
    neighs = copy.(fill(Int[], n))
    weights = copy.(fill(Float64[], n))
    nneighs = zeros(Int,n)

    if criterion == :Queen critthr = 1 end
    if criterion == :Rook  critthr = 2 end

    # Use multithreading to calculate polygons hits
    neighssingle = copy.(fill(Int[], n))
    Threads.@threads for i in 1:n
        @inbounds for j in candidates[i]
            p1 = [x[i], y[i]]
            p2 = [x[j], y[j]]
            polhits = _hits(p1, p2, critthr, tol)
            
            if polhits  
                push!(neighssingle[i], j) 
            end
        end        
    end
    
    # Makes neighbors bilateral
    for i in 1:n
        for j in neighssingle[i]
            push!(neighs[i], j)
            push!(neighs[j], i)

            nneighs[i] += 1
            nneighs[j] += 1
        end
    end

    # Sort neighbours lists and compute weights
    for i in 1:n
        sort!(neighs[i])
        weights[i] = ones(nneighs[i]) ./ nneighs[i]
    end

    SpatialWeights(n, neighs, weights, nneighs, :row)
end

"""
    polyneigh(A, criterion = :Queen)
Build a spatial weights object from table A that contains a geometry column.

# Optional Arguments
- `criterion=:Queen`: neighbour criterion. `:Queen` or `:Rook`.
- `tol=0.0`: tolerance for polygon contiguity.
"""
function polyneigh(A::Any; criterion::Symbol = :Queen, tol::Float64 = 0.0)::SpatialWeights
    
    istable(A) || throw(ArgumentError("Unknown geometry or not polygon geometry"))

    (:geometry in propertynames(A)) || throw(ArgumentError("table does not have :geometry information"))

    return polyneigh(A.geometry, criterion = criterion, tol = tol)
end
