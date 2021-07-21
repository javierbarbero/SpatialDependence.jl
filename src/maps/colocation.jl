# Co-location maps classificator
struct CoLocation <:  AbstractCoLocationMapClassificator
    mcr::AbstractMapClassificator
    CoLocation(mcr::AbstractMapClassificator) = return new(mcr);
end

# Structure for Co-location  Map Classification
struct CoLocationMapClassification <: AbstractMapClassification
    mcr::AbstractMapClassificator
    k::Int
    group::Vector{Int}
    ngroup::Vector{Int}
    grouplabs::Vector{String}
end

function mapclassify(mcr::AbstractCoLocationMapClassificator, x::Vector{Vector{T}} where T;
    lower::Symbol = :open, upper::Symbol = :closed)::CoLocationMapClassification  

    l = length(x)
    n = length(x[1])
    lgroup = zeros(Int, n, l)
    lk = zeros(Int, l)
    comcr = mcr.mcr

    if isa(comcr, AbstractGraduatedMapClassificator)

        # Get categories for each variable
        for i in 1:l
            ni = length(x[i])
            (n == ni) || throw(DimensionMismatch("dimensions must match: x[1] has ($(n)), x[$(i)] has ($ni)"))

            mc = mapclassify(comcr, x[i], lower = lower, upper = upper)
            lgroup[:,i] = assignments(mc)
            lk[i] = length(mc)

            (lk[1] == lk[i]) || throw(DimensionMismatch("classifications categories must match: 1 has ($(lk[1])), $(i) has ($(lk[i]))"))
        end

        k = lk[1] + 1

        # Get common categories
        group = zeros(Int, n)
        ngroup = zeros(Int, k)

        for i in 1:n
            gi = lgroup[i,1]
            if all(lgroup[i,:] .== gi)
                group[i] = gi + 1
                ngroup[gi + 1] += 1
            else
                group[i] = 1
                ngroup[1] += 1
            end
        end

        grouplabs = vcat("Other", string.(1:k-1))
    elseif isa(comcr, AbstractUniqueMapClassificator)
        # Get all categories
        cats = unique(mapreduce(unique, vcat, x))
        sort!(cats)
        k = length(cats) + 1
        group = zeros(Int, n)
        ngroup = zeros(Int, k)
        grouplabs = vcat("Other", string.(cats))

        for i in 1:n
            xi = x[1][i]
            for j in 2:l
                # Count number of coincidences
                c = 1
                if x[j][i] == xi
                    c += 1
                end

                if c == l
                    g = findall(y -> y == xi, cats)[1]
                    group[i] = g + 1
                    ngroup[g + 1] += 1
                else
                    group[i] = 1
                    ngroup[1] += 1
                end
            end
        end

    end

    return CoLocationMapClassification(mcr, k, group, ngroup, grouplabs)
end

function maplabels(mc::CoLocationMapClassification; counts::Bool = true)::Vector{String}

    ngroup = mc.ngroup
    grouplabs = mc.grouplabs
    labels = grouplabs
    labels[1] = "Other"

    if counts
        labels = labels .* " (" .* string.(ngroup) .* ")"
    end

    labels
end

levels(mc::CoLocationMapClassification)::Vector{String} = return mc.grouplabs;
