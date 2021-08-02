# This file contains functions for creating lattices

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
        geometries = MultiPolygon[]
    elseif gtype == :point
        geometries = Point[]
    end

    for i in iiter
        for j in jiter  
            if gtype == :polygon          
                push!(geometries, 
                    MultiPolygon(
                        Polygon([
                            [[0.0 + j, 0.0 + i], 
                            [1.0 + j, 0.0 + i], 
                            [1.0 + j, 1.0 + i], 
                            [0.0 + j, 1.0 + i] ]])))
            elseif gtype == :point
                push!(geometries, Point([0.0 + j, 0.0 + i]))
            end
        end
    end

    return geometries
end