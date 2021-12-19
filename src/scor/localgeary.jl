# Local Geary test of Spatial Autocorrelation
struct LocalGeary <: AbstractLocalSpatialAutocorrelation
    n::Int
    C::Vector{Float64}
    p::Vector{Float64}
    Cperms::Matrix{Float64}
    Cpermsmean::Vector{Float64}
    Cpermsstd::Vector{Float64}
    z::Vector{Float64}
    q::Vector{Symbol}
    categories::Symbol
end

"""
    localgeary(x, W)
Compute the Local Geary test of spatial autocorrelation.

# Optional Arguments
- `permutations=9999`: number of permutations for the randomization test.
- `rng=default_rng()`: random number generator for the randomization test.
- `corrected=true`: divide the scaling factor by ``n-1`` instead of ``n``.
- `categories=:positivenegative`: assing observations to positive or negative spatial autocorrelation, or in combination with the `:moran` scatterplot.
"""
function localgeary(x::AbstractVector{T} where T<:Real, W::SpatialWeights; permutations::Int = 9999,
    corrected::Bool = true, categories::Symbol = :positivenegative,
    rng::AbstractRNG = default_rng())::LocalGeary

    (categories == :positivenegative) || (categories == :moran) || throw(ArgumentError("`categories` must be :positivenegative or :moran"))

    n = length(x)
    z = x .- mean(x)   
    
    m2::Float64 = sum(z.^2) 
    if corrected
        m2 = m2 ./ (n - 1)
    else
        m2 = m2 ./ n
    end

    # Auxiliar function to calculate Local Moran
    function localgeary_calc(zi::Number, wi::AbstractVector, zneighi::AbstractVector)::Float64
        return (1 ./ m2)  .* sum(wi .* (zi .- zneighi).^2 )
    end
    
    # Local Geary
    C = zeros(n)
    for i in 1:n
        zi = z[i]
        wi = weights(W, i)
        zneighi = z[neighbors(W, i)]
        C[i] = localgeary_calc(zi, wi, zneighi) 
    end

    # Conditional randomizatoin
    Cperms = crand_local(permutations, z, W, localgeary_calc, rng)
    
    larger = sum(Cperms .>= repeat(C, 1, permutations), dims = 2)
    low = (permutations .- larger) .< larger
    larger[low] .= permutations .- larger[low]
    p = (larger .+ 1) ./ (permutations + 1)
    p = vec(p)

    
    Cpermsstd = vec(std(Cperms, dims = 2, corrected = false))
    Cpermsmean = vec(mean(Cperms, dims = 2))
    zval = (C .- Cpermsmean)  ./ Cpermsstd
    zval = vec(zval)

    # Classification
    Cmean = mean(C)
    Wz = slag(W, z)

    q = Array{Symbol}(undef, n)
    if categories == :positivenegative
        for i in 1:n
            if zval[i] < 0
                q[i] = :P
            else
                q[i] = :N
            end
        end
    elseif categories == :moran
        for i in 1:n
            if zval[i] < 0
                if (z[i] > 0) & (Wz[i] > 0)
                    q[i] = :HH
                elseif (z[i] < 0) & (Wz[i] < 0)
                    q[i] = :LL
                else
                    q[i] = :OP
                end
            else
                q[i] = :NE
            end
        end
    end

    return LocalGeary(n, C, p, Cperms, Cpermsmean, Cpermsstd, zval, q, categories)

end

score(x::LocalGeary) = x.C;

scoreperms(x::LocalGeary) = x.Cperms;

mean(x::LocalGeary) = x.Cpermsmean;

std(x::LocalGeary) = x.Cpermsstd;

zscore(x::LocalGeary) = x.z;

pvalue(x::LocalGeary) = x.p;

assignments(x::LocalGeary) = x.q;

testname(::LocalGeary) = "Local Geary";

function labelsorder(x::LocalGeary) 
    if x.categories == :positivenegative
        return [:ns, :P, :N]
    elseif x.categories == :moran
        return [:ns; :HH; :LL; :OP; :NE];
    end
end

function labelsnames(x::LocalGeary) 
    if x.categories == :positivenegative
        return ["Not Significant"; "Positive"; "Negative"]
    elseif x.categories == :moran
        return ["Not Significant"; "High-High"; "Low-Low"; "Other Positive"; "Negative"]
    end
end

function labelcolor(::LocalGeary, x::String)

    # :positivenegative categories
    if x == "P"
        catcolor = :red
        alpha = 1.0
    elseif x == "N"
        catcolor =:blue
        alpha = 1.0
    # :moran categories
    elseif x == "HH"
        catcolor = "#B2182B" # (178, 24, 43)
        alpha = 1.0
    elseif x == "LL"
        catcolor = "#EF8A62" # (239, 138, 98)
        alpha = 1.0
    elseif x == "OP"
        catcolor = "#FDDBC7" # (253, 219, 199)
        alpha = 1.0
    elseif x == "NE"
        catcolor = "#67ADC7" # (103, 173, 199)
        alpha = 1.0
    # not significant
    elseif x == "ns"
        catcolor = :lightgrey
        alpha = 0.4
    end

    return (catcolor, alpha)
end