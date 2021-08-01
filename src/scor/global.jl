# Abstract type for Global Spatial Autocorrelation
abstract type AbstractGlobalSpatialAutocorrelation end

function Base.show(io::IO, x::AbstractGlobalSpatialAutocorrelation)

    println(io, testname(x) * " test of Global Spatial Autocorrelation")
    println(io, "--------------------------------------------")
    println(io, "")
    println(io, testname(x) * ": ", round(score(x), digits = 7))
    println(io, "Expectation: ", round(expected(x), digits = 7))
    println(io, "")
    println(io, "Randomization test with ", length(scoreperms(x)), " permutations.")
    println(io, " Mean: ", round(mean(x), digits = 7))
    println(io, " Std Error: ", round(std(x), digits = 7))
    println(io, " zscore: ", round(zscore(x), digits = 7))
    println(io, " p-value: ", round(pvalue(x), digits = 7))

end
