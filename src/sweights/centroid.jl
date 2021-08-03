# This file contains functions to calculate polygon centroid and mean center

"""
    meancenter(P)
Calculate the mean center of a polygon or vector of polygons ``P``.
"""
function meancenter end

function meancenter(P::Union{AbstractPolygon,AbstractMultiPolygon})::Tuple{Float64, Float64}
    xy = _getpointsPoligon(P)

    cx = mean(xy[1])
    cy = mean(xy[2])

    return cx, cy
end

function meancenter(P::Vector{T} where T <:Union{Missing,AbstractPolygon,AbstractMultiPolygon})::Tuple{Vector{Float64}, Vector{Float64}}
    cxy = meancenter.(P)

    cx = first.(cxy)
    cy = last.(cxy)

    return cx, cy
end

"""
    centroid(P)
Calculate the centroid of a polygon or vector of polygons ``P``.
"""
function centroid end

function centroid(P::Union{AbstractPolygon,AbstractMultiPolygon})::Tuple{Float64, Float64}
    xy = _getpointsPoligon(P)

    xo = xy[1]
    yo = xy[2]

    # Formula: https://en.wikipedia.org/wiki/Centroid#Of_a_polygon

    # If the last point is not equal to the first point, add a new point similar to the first point.
    if (xo[end] != xo[1]) | (yo[end] != yo[1])
        x = copy(xo)
        y = copy(yo)
        
        push!(x, x[1])
        push!(y, y[1])    
    else
        x = xo
        y = yo
    end

    # Loop to calculate elements of A and centroids
    n = length(x)
    A = 0
    xprodsum = 0
    yprodsum = 0

    for i in 1:n-1
        xycross = x[i] * y[i + 1] - x[i + 1] * y[i]
        xprod = (x[i] + x[i + 1]) * xycross
        yprod = (y[i] + y[i + 1]) * xycross

        A += xycross
        xprodsum += xprod
        yprodsum += yprod
    end
    A = A / 2

    # Centroids
    cx = 1 / (6 * A) * xprodsum
    cy = 1 / (6 * A) * yprodsum

    return cx, cy
end

function centroid(P::Vector{T} where T <:Union{Missing,AbstractPolygon,AbstractMultiPolygon})::Tuple{Vector{Float64}, Vector{Float64}}
    cxy = centroid.(P)

    cx = first.(cxy)
    cy = last.(cxy)

    return cx, cy
end
