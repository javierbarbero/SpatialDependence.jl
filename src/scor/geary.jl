# Geary's c test of Global Spatial Autocorrelation
struct GlobalGeary <: AbstractGlobalSpatialAutocorrelation
    n::Int
    C::Float64
    EC::Float64
    p::Float64
    Cperms::Vector{Float64}
    Cpermsmean::Float64
    Cpermsstd::Float64
    z::Float64
end

"""
    geary(x, W)
Compute the Geary's c test of global spatial autocorrelation.

# Optional Arguments
- `permutations=9999`: number of permutations for the randomization test.
- `rng=default_rng()`: random number generator for the randomization test.
"""
function geary(x::AbstractVector{T} where T<:Real, W::SpatialWeights; permutations::Int = 9999,
    rng::AbstractRNG = default_rng())::GlobalGeary

    n = length(x)
    z = x    
    S0 = sum(sum.(W.weights::Vector{Vector{Float64}}))
    denominator = sum( (z .- mean(z)).^2 )

    # Auxiliar function to calculate Geary's C
    function geary_calc(z::AbstractVector, W::SpatialWeights)::Float64
        numerator = 0.0       
        for i in 1:n
            numerator += sum( weights(W, i) .* (z[i] .- z[neighbors(W, i)]).^2 )
        end
        return (n - 1)/(2 * S0) * numerator / denominator
    end
    
    # Gearys's C
    C = geary_calc(z, W)
    EC = 1.0

    # Build permutations array
    Zi = Array{Float64}(undef, n, permutations)
    for i in 1:permutations
        Zi[:,i] = shuffle(rng, z)
    end

    # Calculate Geary's c for all the permutations
    Cperms = zeros(permutations)
    Threads.@threads for i in 1:permutations
        Cperms[i] = geary_calc(view(Zi, :, i), W)
    end

    larger = sum(Cperms .>= C)
    if (permutations - larger) < larger 
        larger = permutations - larger
    end
    p = (larger + 1) / (permutations + 1)
    
    Cpermsstd = std(Cperms, corrected = false)
    Cpermsmean = mean(Cperms)
    z = (C - Cpermsmean)  / Cpermsstd

    return GlobalGeary(n, C, EC, p, Cperms, Cpermsmean, Cpermsstd, z)

end

score(x::GlobalGeary) = x.C;

scoreperms(x::GlobalGeary) = x.Cperms;

expected(x::GlobalGeary) = x.EC;

mean(x::GlobalGeary) = x.Cpermsmean;

std(x::GlobalGeary) = x.Cpermsstd;

zscore(x::GlobalGeary) = x.z;

pvalue(x::GlobalGeary) = x.p;

testname(x::GlobalGeary) = "Geary's c";
