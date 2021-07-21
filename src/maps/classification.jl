# This file contains types and functions for Choropleth Maps

# Abstract Map Classificator
abstract type AbstractMapClassificator end

abstract type AbstractGraduatedMapClassificator <: AbstractMapClassificator end

abstract type AbstractStatisticalMapClassificator  <: AbstractGraduatedMapClassificator end

abstract type AbstractUniqueMapClassificator <: AbstractMapClassificator end

abstract type AbstractCoLocationMapClassificator <: AbstractMapClassificator end

# Default color palettes
# This functions returns the default color palette and a boolean indicating if colors are reversed 
defaultpalette(::AbstractGraduatedMapClassificator) = (:YlOrBr, false)

defaultpalette(::AbstractStatisticalMapClassificator) = (:RdBu, true)

defaultpalette(::AbstractUniqueMapClassificator) = (:Paired, false)

defaultpalette(mcr::AbstractCoLocationMapClassificator) = defaultpalette(mcr.mcr)

# Abstract Map Classification
abstract type AbstractMapClassification end

"""
    mapclassify (mcr::AbstractMapClassificator, x::Vector)
Classify observations in variable `x` in classess using the criterion specified in `mcr`.
"""
function mapclassify end

"""
    maplabels (mcAbstractMapClassification)
Get the labels of the classification `mc`.
# Optional Arguments
- `digits=2`: number of decimal digits.
- `sep=", "`: separator between class lower and upper bounds.
- `counts=true`: include the total number of observations on each class.
"""
function maplabels end

"""
    counts(mc::AbstractMapClassification)
Get the vector of classess sizes.

`counts(mc)[k]` is the number of observations assigned to the ``k``-class.
"""
counts(mc::AbstractMapClassification)::Vector{Int} = return mc.ngroup;

"""
    bounds(mc::AbstractMapClassification)
Get the lower and upper bounds of the classes.

`bounds(mc)[1][k]` is the lower bound of the ``k``-class.

`bounds(mc)[2][k]` is the upper bound of the ``k``-class.
"""
function bounds(::AbstractMapClassification)::Tuble{Vector{Float64}, Vector{Float64}} end

"""
    levels(mc::AbstractMapClassification)
Get the levels of the classes.

`levels(mc)[k]` is the level of the ``k``-class.
"""
function levels(::AbstractMapClassification)::Vector{String} end

"""
    assignments(mc::AbstractMapClassification)
Get the vector of classes indices for each observation.

`assignments(mc)[i]` is the index of the class to which the ``i``-th observation is assigned.
"""
assignments(mc::AbstractMapClassification)::Vector{Int} = return mc.group;

length(mc::AbstractMapClassification)::Int = return mc.k;
