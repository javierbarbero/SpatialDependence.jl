# Map classificator for Unique values
struct Unique <: AbstractUniqueMapClassificator 
    uniq::Union{Vector{String}, Vector{Symbol}, Nothing}
    Unique(uniq::Union{Vector{String}, Vector{Symbol}, Nothing} = nothing) = return new(uniq);
end

# Structure for Unique Values  Map Classification
struct UniqueMapClassification <: AbstractMapClassification
    mcr::AbstractMapClassificator
    k::Int
    group::Vector{Int}
    ngroup::Vector{Int}
    grouplabs::Vector{String}
end

function mapclassify(mcr::AbstractUniqueMapClassificator, x::AbstractVector)::UniqueMapClassification  

    if isnothing(mcr.uniq)
        uniq = unique(x)
    else
        uniq = mcr.uniq
    end
    k = length(uniq)
    group = zeros(Int, length(x))
    ngroup = zeros(Int, k)
    grouplabs = string.(uniq)

    for i in 1:k
        subset = x .== uniq[i]
        group[subset] .= i
        ngroup[i] = sum(subset)
    end

    return UniqueMapClassification(mcr, k, group, ngroup, grouplabs)
end

function maplabels(mc::UniqueMapClassification; counts::Bool = true)::Vector{String}

    ngroup = mc.ngroup
    grouplabs = mc.grouplabs
    labels = grouplabs

    if counts
        labels = labels .* " (" .* string.(ngroup) .* ")"
    end

    labels
end

levels(mc::UniqueMapClassification)::Vector{String} = return mc.grouplabs;
