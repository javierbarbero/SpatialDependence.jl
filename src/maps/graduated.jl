# Common graduated values classification
struct EqualIntervals <:  AbstractGraduatedMapClassificator 
    k::Int
    EqualIntervals(k::Int = 5) = return new(k);
end

struct Quantiles <: AbstractGraduatedMapClassificator 
    k::Int
    Quantiles(k::Int = 5) = return new(k);
end

struct NaturalBreaks <: AbstractGraduatedMapClassificator 
    k::Int
    NaturalBreaks(k::Int = 5) = return new(k);
end

struct CustomBreaks <: AbstractGraduatedMapClassificator 
    k::Int
    bins::Union{Vector{Int}, Vector{Float64}}
    CustomBreaks(bins::Union{Vector{Int}, Vector{Float64}}) = return new(length(bins) + 1, bins);
end

# Bounds
function mapclassifybounds(mcr::EqualIntervals, x::Vector{T} where T)::Tuple{Vector{Float64}, Vector{Float64}}
    k = mcr.k 
    minx = minimum(x)
    maxx = maximum(x)
    stepinc = (maxx - minx) / k 

    lbound = zeros(k)
    ubound = zeros(k)

    for i in 1:k
        lbound[i] = minx + (i - 1) * stepinc
        ubound[i] = minx + i * stepinc
    end

    lbound, ubound
end

function mapclassifybounds(mcr::Quantiles, x::Vector{T} where T)::Tuple{Vector{Float64}, Vector{Float64}}
    k = mcr.k
    qx = quantile(x, range(1/k, 1, length = k))

    lbound = zeros(k)
    lbound[1] = minimum(x)
    lbound[2:k] = qx[1:k-1]

    ubound = qx

    lbound, ubound
end

function mapclassifybounds(mcr::CustomBreaks, x::Vector{T} where T)::Tuple{Vector{Float64}, Vector{Float64}}
    bins = mcr.bins

    lbound = vcat(minimum(x), bins) 
    ubound = vcat(bins, maximum(x))

    lbound, ubound
end


function mapclassifybounds(mcr::NaturalBreaks, x::Vector{T} where T)::Tuple{Vector{Float64}, Vector{Float64}}
    k = mcr.k
    n = length(x)
    xs = sort(x)

    (k < n) || throw(ArgumentError("number of classes should be lower than number of observations"))

    # Create matrices and set initial values
    mat1 = zeros(Int64,   n, k) # Lower class limits
    mat2 = zeros(Float64, n, k) # Variance combinations

    mat1[1, :]     .= 1
    mat2[2:end, :] .= Inf   # Set to the most larger possible value 

    # Loop through each observation
    v = 0.0
    @inbounds for l in 2:n
        s1 = 0.0 # sum
        s2 = 0.0 # sum of squares
        w = 0 # Number of observations considered

        # Loop through observations till l in reverse order
        for m in reverse(1:l)
            val = xs[m]

            # Increase sum of sqaures and observations 
            s1 = s1 + val
            s2 = s2 + val * val
            w = w + 1

            # Variance
            v = s2 - (s1 * s1) / w

            # Update lower bound and variance if newer is lower
            prev = m - 1
            if prev != 0
                for j in 2:k
                    if mat2[l, j] >= (v + mat2[prev, j - 1])                       
                        mat1[l, j] = m
                        mat2[l, j] = v + mat2[prev, j - 1]
                    end
                end
            end
        end

        # Store lower class limits and variance
        mat1[l, 1] = 1
        mat2[l, 1] = v

    end

    # Get class bounds
    kclass = zeros(k + 1)
    kclass[end] = xs[end]
    kclass[1] = xs[1]

    # Get class bounds
    lbound = zeros(k)
    ubound = zeros(k)
    
    pivot = n
    for ki in reverse(2:k)
        pivot = mat1[pivot, ki] - 1
        lbound[ki] = xs[pivot]
    end
    
    lbound[1] = minimum(xs)
    ubound[1:k-1] = lbound[2:k]
    ubound[k] = maximum(xs)

    lbound, ubound
end

# Structure for Graduated  Map Classification
struct GraduatedMapClassification <: AbstractMapClassification
    mcr::AbstractMapClassificator
    k::Int
    group::Vector{Int}
    ngroup::Vector{Int}
    lbound::Vector{Float64}
    ubound::Vector{Float64}
    lower::Symbol
    upper::Symbol
end

function mapclassify(mcr::AbstractGraduatedMapClassificator, x::Vector{T} where T; 
    lower::Symbol = :open, upper::Symbol = :closed)::GraduatedMapClassification  

    ((lower == :open) || (lower == :closed )) || throw(ArgumentError("lower bound must be :open or :closed"))
    ((upper == :open) || (upper == :closed )) || throw(ArgumentError("upper bound must be :open or :closed"))
    (lower !== upper) || throw(ArgumentError("lower and upper bounds cannot be both :open or :closed"))

    lbound, ubound = mapclassifybounds(mcr, x)

    k = mcr.k
    group = zeros(Int, length(x))
    ngroup = zeros(Int, k)

    for i in 1:k        
        if i == 1
            # First group
            if upper == :open
                subset = (x .>= lbound[i]) .& (x .< ubound[i])
            elseif upper == :closed
                subset = (x .>= lbound[i]) .& (x .<= ubound[i])
            end
        elseif i == k
            # Last group
            if lower == :open
                subset = (x .> lbound[i]) .& (x .<= ubound[i])
            else
                subset = (x .>= lbound[i]) .& (x .<= ubound[i])
            end            
        else i < k
            # Intermediate groups
            if lower == :open && upper == :closed
                subset = (x .> lbound[i]) .& (x .<= ubound[i])
            elseif lower == :closed && upper == :open
                subset = (x .>= lbound[i]) .& (x .< ubound[i])
            end
        end

        group[subset] .= i
        ngroup[i] = sum(subset)
    end

    return GraduatedMapClassification(mcr, k, group, ngroup, lbound, ubound, lower, upper)
end

function maplabels(mc::GraduatedMapClassification; digits::Int = 2, sep::String = ", ", counts::Bool = true)::Vector{String}
    k = mc.k
    ngroup = mc.ngroup
    lbound = mc.lbound
    ubound = mc.ubound
    lower = mc.lower
    upper = mc.upper

    labels = Vector{String}(undef, k)

    # Label for each group in the classification
    for i in 1:k
        labeli = ""

        # Lower interval Symbol
        if i == 1
            labeli = labeli * "["
        else
            if lower == :open
                labeli = labeli * "("
            elseif lower == :closed
                labeli = labeli * "["
            end
        end
 
        # Range
        labeli = labeli * string(round(lbound[i], digits = digits)) * 
            sep * 
            string(round(ubound[i], digits = digits))

        # Upper interval symbol
        if i == k
            labeli = labeli * "]"
        else
            if upper == :open
                labeli = labeli * ")"
            elseif upper == :closed
                labeli = labeli * "]"
            end
        end

        # Counts
        if counts
            labeli = labeli * " (" * string(ngroup[i]) * ")"
        end

        labels[i] = labeli
    end

    labels
end

bounds(mc::GraduatedMapClassification)::Tuple{Vector{Float64}, Vector{Float64}} = return mc.lbound, mc.ubound;

levels(mc::GraduatedMapClassification)::Vector{String} = return string.(1:mc.k);
