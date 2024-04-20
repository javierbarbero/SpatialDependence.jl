# Getis-Ord statistic
struct GetisOrd <: AbstractLocalSpatialAutocorrelation
    n::Int
    G::Vector{Float64}
    p::Vector{Float64}
    Gperms::Matrix{Float64}
    Gpermsmean::Vector{Float64}
    Gpermsstd::Vector{Float64}
    z::Vector{Float64}
    q::Vector{Symbol}
    star::Bool
end

"""
    getisord(x, W)
Compute the Getis-Ord statistic.

# Optional Arguments
- `star=true`: compute the Gi* statistic, or the Gi if set to `false`.
- `permutations=9999`: number of permutations for the randomization test.
- `rng=default_rng()`: random number generator for the randomization test.
"""
function getisord(x::AbstractVector{T} where T, W::SpatialWeights; permutations::Int = 9999,
    star::Bool = true, rng::AbstractRNG = default_rng())::GetisOrd

    wt = wtransformation(W) 
    wt == :row || wt == :binary || throw(ArgumentError("W must be row standardized or binary"))

    n = length(x)
    
    denon::Float64 = sum(x)

    # Auxiliar function to calculate Getis Ord
    function getisord_calc(xi::Number, wi::AbstractVector, xneighi::AbstractVector)::Float64
        return sum( wi .* xneighi) ./ (denon - xi)
    end

    function getisord_calc_star(xi::Number, wi::AbstractVector, xneighi::AbstractVector)::Float64
        if wt == :row
            wistar = 1 / (length(wi) + 1)
        elseif wt == :binary
            wistar = 1
        end
        return (wistar .* xi + sum(wistar .* xneighi)) ./ (denon)
    end

    if star
        getisord_calc_fun = getisord_calc_star
    else
        getisord_calc_fun = getisord_calc
    end
    
    # Getis-Ord
    G = zeros(n)
    for i in 1:n
        xi = x[i]
        wi = weights(W, i)
        xneighi = x[neighbors(W, i)]
        G[i] = getisord_calc_fun(xi, wi, xneighi) 
    end

    # Conditional randomizatoin
    Gperms = crand_local(permutations, x, W, getisord_calc_fun, rng)
    
    larger = sum(Gperms .>= repeat(G, 1, permutations), dims = 2)
    low = (permutations .- larger) .< larger
    larger[low] .= permutations .- larger[low]
    p = (larger .+ 1) ./ (permutations + 1)
    p = vec(p)

    
    Gpermsstd = vec(std(Gperms, dims = 2, corrected = false))
    Gpermsmean = vec(mean(Gperms, dims = 2))
    zval = (G .- Gpermsmean)  ./ Gpermsstd
    zval = vec(zval)

    # Classification
    q = Array{Symbol}(undef, n)
    for i in 1:n
        if zval[i] > 0
            q[i] = :H
        else
            q[i] = :L
        end
    end

    return GetisOrd(n, G, p, Gperms, Gpermsmean, Gpermsstd, zval, q, star)

end

score(x::GetisOrd) = x.G;

scoreperms(x::GetisOrd) = x.Gperms;

mean(x::GetisOrd) = x.Gpermsmean;

std(x::GetisOrd) = x.Gpermsstd;

zscore(x::GetisOrd) = x.z;

pvalue(x::GetisOrd) = x.p;

assignments(x::GetisOrd) = x.q;

testname(x::GetisOrd) = x.star ? "Getis-Ord Gi*" : "Getis-Ord Gi";

function labelsorder(::GetisOrd) 
    return [:ns, :H, :L]
end

function labelsnames(::GetisOrd) 
    return ["Not Significant"; "High"; "Low"]
end

function labelcolor(::GetisOrd, x::String)

    if x == "H"
        catcolor = :red
        alpha = 1.0
    elseif x == "L"
        catcolor =:blue
        alpha = 1.0
    # not significant
    elseif x == "ns"
        catcolor = :lightgrey
        alpha = 0.4
    end

    return (catcolor, alpha)
end