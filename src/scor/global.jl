# Abstract type for Global Spatial Autocorrelation
abstract type AbstractGlobalSpatialAutocorrelation end

function Base.show(io::IO, x::AbstractGlobalSpatialAutocorrelation)

    println(io, testname(x) * " test of Global Spatial Autocorrelation")
    println(io, "--------------------------------------------")
    println(io, "")
    println(io, testname(x) * ": ", round(score(x), digits = 7))
    println(io, "Expectation: ", round(expected(x), digits = 7))
    println(io, "")

    permutations = length(scoreperms(x))
    println(io, "Randomization test with ", permutations, " permutations.")

    if permutations > 0
        println(io, " Mean: ", round(mean(x), digits = 7))
        println(io, " Std Error: ", round(std(x), digits = 7))
        println(io, " zscore: ", round(zscore(x), digits = 7))
        println(io, " p-value: ", round(pvalue(x), digits = 7))
    else
        println(io, "Statistical significance cannot be assessed with 0 permutataions.")
    end

end
