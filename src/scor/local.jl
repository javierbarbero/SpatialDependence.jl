# Abstract type for Local Spatial Autocorrelation
abstract type AbstractLocalSpatialAutocorrelation end

"""
    issignificant(x, α, adjust = :none)
Return a vector of boolean values indicating if the local statistics are significant or not at the desired threshold ``α``.

p-values can be adjusted with the `adjust` parameter using the Bonferroni correction `:bonferroni` or controlling for the False Discovery Rate, `:fdr`
"""
function issignificant(x::AbstractLocalSpatialAutocorrelation, α::Float64; adjust::Symbol = :none)::Vector{Bool}
    p = pvalue(x)
    n = length(p)

    (adjust == :none || adjust == :bonferroni || adjust == :fdr) ||  
        throw(ArgumentError("unknown p-value adjustment $(adjust)"))

    if adjust == :none
        return p .< α
    elseif adjust == :bonferroni
        return p .< (α / n)
    elseif adjust == :fdr
        psort = sort(p)
        pfdr = (1:n) .* α ./ n
        lower = psort .< pfdr
        threshold = findfirst(x -> x == false, lower)
        return p .< pfdr[threshold]
    end
end

function Base.show(io::IO, x::AbstractLocalSpatialAutocorrelation)

    println(io, testname(x) * " test of Spatial Autocorrelation")
    println(io, "--------------------------------------------")
    println(io, "")
    println(io, "Randomization test with ", size(scoreperms(x), 2), " permutations.")
    println(io, "Interesting locations at 0.05 significance level:")

    # Display clusters
    q = assignments(x)
    p = pvalue(x)

    labels = labelsorder(x)
    nl = length(labels)
    countcat = zeros(Int, nl)
    for i in 1:nl
        issig = issignificant(x, 0.05, adjust = :fdr)
        if labels[i] == :ns
            countcat[i] = count(.! issig)
        else
            countcat[i] = count(issig .& (q .== labels[i]))
        end
    end

    labelsstr = labelsnames(x)
    labelmaxlength = maximum(length.(labelsstr))
    for i in 1:nl
        if labels[i] != :ns
            println(io, " ", lpad(labelsstr[i], labelmaxlength), ": ", countcat[i])
        end
    end

end
