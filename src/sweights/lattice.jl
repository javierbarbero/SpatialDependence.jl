# This file contains functions for creating lattices

# GeoInterface Point Implementation
struct LatticePoint
    x::Float64
    y::Float64
end

function LatticePoint(v::Vector{Float64})
    return LatticePoint(v[1], v[2])
end

GI.isgeometry(::LatticePoint) = true
GI.geomtrait(::LatticePoint) = GI.PointTrait()
GI.ncoord(::GI.PointTrait, geom::LatticePoint) = 2
GI.getcoord(::GI.PointTrait, geom::LatticePoint, i) = [geom.x, geom.y][i]

# GeoInterface LineStringTrait Implementation
struct LatticeCurve
    p::Vector{LatticePoint}
end

GeoInterface.isgeometry(::LatticeCurve) = true
GeoInterface.geomtrait(::LatticeCurve) = GI.LineStringTrait()
GeoInterface.ngeom(::GI.LineStringTrait, geom::LatticeCurve) = length(geom.p)
GeoInterface.getgeom(::GI.LineStringTrait, geom::LatticeCurve, i) = geom.p[i]
Base.convert(T::Type{LatticeCurve}, geom::X) where {X} = Base.convert(T, GI.geomtrait(geom), geom)
Base.convert(::Type{LatticeCurve}, ::GI.LineStringTrait, geom::LatticeCurve) = geom

function LatticeCurve(c::Vector{Vector{Float64}})
    return LatticeCurve(LatticePoint.(c))
end

# GeoInterface Polygon Implementation
struct LatticePolygon
    c::LatticeCurve
end

function LatticePolygon(c::Vector{Vector{Float64}})
    return LatticePolygon(LatticeCurve(c))
end

GI.isgeometry(::LatticePolygon) = true
GI.geomtrait(::LatticePolygon) = GI.PolygonTrait()
GI.ngeom(::GI.PolygonTrait, geom::LatticePolygon) = 1
GI.getgeom(::GI.PolygonTrait, geom::LatticePolygon, i) = geom.c

GI.is3d(::GI.PolygonTrait, ::LatticePolygon) = false

# GeoInterface MultiPolygon Implementation
struct LatticeMultiPolygon
    p::Vector{LatticePolygon}
end

function LatticeMultiPolygon(p::LatticePolygon)
    return LatticeMultiPolygon([p])
end

function LatticeMultiPolygon(p::Vector{Vector{Vector{Vector{Float64}}}})
    return LatticeMultiPolygon([LatticePolygon.(p)])
end

GI.isgeometry(::LatticeMultiPolygon) = true
GI.geomtrait(::LatticeMultiPolygon) = GI.MultiPolygonTrait()
GI.ngeom(::GI.MultiPolygonTrait, geom::LatticeMultiPolygon) = length(geom.p)
GI.getgeom(::GI.MultiPolygonTrait, geom::LatticeMultiPolygon, i) = geom.p[i]

GI.is3d(::GI.MultiPolygonTrait, ::LatticeMultiPolygon) = false

"""
    reggeomlattice(n, m)
Create a regular geometry-lattice of ``n`` rows and ``m`` columns.

# Optional Arguments
- `gtype=:polygon`: geometry. :polygon or :point.
- `direction=:rightdown`: direction. :rightup, :rightdown, :leftdown or :leftup.
"""
function reggeomlattice(n::Int, m::Int; gtype::Symbol = :polygon, direction::Symbol = :rightdown)
    n > 0 || throw(ArgumentError("The number of rows must be greater than 0"))
    m > 0 || throw(ArgumentError("The number of columns must be greater than 0"))

    (direction == :rightup) || (direction == :rightdown) || 
    (direction == :leftdown) || (direction == :leftup) ||
        throw(ArgumentError("`direction` must be :rightup, :rightdown, :leftdown or :leftup"))

    (gtype == :polygon) || (gtype == :point) || throw(ArgumentError("`gtype` must be :polygon or :point"))

    if direction == :rightup
        iiter = 1:n
        jiter = 1:m
    elseif direction == :rightdown
        iiter = reverse(1:n)
        jiter = 1:m
    elseif direction == :leftdown
        iiter = reverse(1:n)
        jiter = reverse(1:m)
    elseif direction == :leftup
        iiter = 1:n
        jiter = reverse(1:m)
    end

    if gtype == :polygon
        geometries = LatticeMultiPolygon[]
    elseif gtype == :point
        geometries = LatticePoint[]
    end

    for i in iiter
        for j in jiter  
            if gtype == :polygon          
                push!(geometries, 
                    LatticeMultiPolygon(
                        LatticePolygon([
                            [0.0 + j, 0.0 + i], 
                            [1.0 + j, 0.0 + i], 
                            [1.0 + j, 1.0 + i], 
                            [0.0 + j, 1.0 + i]])))
            elseif gtype == :point
                push!(geometries, LatticePoint(0.0 + j, 0.0 + i))
            end
        end
    end

    return geometries
end