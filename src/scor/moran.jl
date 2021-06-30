
struct GlobalMoran
    n::Int
    I::Float64
    EI::Float64
    p::Float64
    Iperms::Vector{Float64}
    Ipermsstd::Float64
    z::Float64
end

"""
    moran(x, W)
Compute the global moran index of spatial autocorrelation.

# Optional Arguments
- `permutations=9999`: number of permutations for the randomization test.
- `rng=default_rng()`: random number generator for the randomization test.
"""
function moran(x::Vector{Float64}, W::SpatialWeights; permutations::Int64 = 9999,
    rng::AbstractRNG = default_rng())::GlobalMoran

    n = length(x)
    z = x .- mean(x)    
    
    S0 = sum(sum.(W.weights))

    # Auxiliar function to calculate Moran's I
    function moran_calc(z::Vector, W::SpatialWeights)
        Wz = slag(W, z)
        return (sum(Wz .* z) / S0) / (sum(z .* z) / n)
    end
    
    # Morans's I
    I = moran_calc(z, W)
    EI = -1.0 / (n - 1)

    # Randomization
    Iperms = zeros(permutations)
    for i in 1:permutations
        zi = shuffle(rng, z)
        Iperms[i] = moran_calc(zi, W)
    end

    larger = sum(Iperms .>= I)
    p = (larger + 1) / (permutations + 1)
    
    Ipermsstd = std(Iperms, corrected = false)
    EIperms = sum(Iperms) / permutations
    z = (I - EIperms)  / Ipermsstd

    return GlobalMoran(n, I, EI, p, Iperms, Ipermsstd, z)

end

score(x::GlobalMoran) = x.I;

std(x::GlobalMoran) = x.Ipermsstd;

zscore(x::GlobalMoran) = x.z;

pvalue(x::GlobalMoran) = x.p;

function Base.show(io::IO, x::GlobalMoran)

    println(io, "Global Moran test of Spatial Autocorrelation")
    println(io, "--------------------------------------------")
    println(io, "")
    println(io, "Moran's I: ", round(score(x), digits = 7))
    println(io, "Expectation:", round(x.EI, digits = 7))
    println(io, "")
    println(io, "Randomization test with ", length(x.Iperms), " permutations.")
    println(io, "Standard Error: ", round(std(x), digits = 7))
    println(io, "zscore: ", round(zscore(x), digits = 7))
    println(io, "p-value: ", round(pvalue(x), digits = 7))

end
