# Moran's I test of Global Spatial Autocorrelation
struct GlobalMoran <: AbstractGlobalSpatialAutocorrelation
    n::Int
    I::Float64
    EI::Float64
    p::Float64
    Iperms::Vector{Float64}
    Ipermsmean::Float64
    Ipermsstd::Float64
    z::Float64
end

"""
    moran(x, W)
Compute the  Moran's I test of global spatial autocorrelation.

# Optional Arguments
- `permutations=9999`: number of permutations for the randomization test.
- `rng=default_rng()`: random number generator for the randomization test.
"""
function moran(x::Vector{T} where T<:Real, W::SpatialWeights; permutations::Int = 9999,
    rng::AbstractRNG = default_rng())::GlobalMoran

    n = length(x)
    z = x .- mean(x)    
    
    S0 = sum(sum.(W.weights::Vector{Vector{Float64}}))

    # Auxiliar function to calculate Moran's I
    function moran_calc(z::Vector, W::SpatialWeights)::Float64
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
    Ipermsmean = mean(Iperms)
    z = (I - Ipermsmean)  / Ipermsstd

    return GlobalMoran(n, I, EI, p, Iperms, Ipermsmean, Ipermsstd, z)

end

score(x::GlobalMoran) = x.I;

scoreperms(x::GlobalMoran) = x.Iperms;

mean(x::GlobalMoran) = x.Ipermsmean;

expected(x::GlobalMoran) = x.EI;

std(x::GlobalMoran) = x.Ipermsstd;

zscore(x::GlobalMoran) = x.z;

pvalue(x::GlobalMoran) = x.p;

testname(x::GlobalMoran) = "Moran's I";
