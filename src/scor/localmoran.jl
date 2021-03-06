# Local Moran test of Spatial Autocorrelation
struct LocalMoran <: AbstractLocalSpatialAutocorrelation
    n::Int
    I::Vector{Float64}
    p::Vector{Float64}
    Iperms::Matrix{Float64}
    Ipermsmean::Vector{Float64}
    Ipermsstd::Vector{Float64}
    z::Vector{Float64}
    q::Vector{Symbol}
end

"""
    localmoran(x, W)
Compute the Local Moran test of spatial autocorrelation.

# Optional Arguments
- `permutations=9999`: number of permutations for the randomization test.
- `rng=default_rng()`: random number generator for the randomization test.
- `corrected=true`: divide the scaling factor by ``n-1`` instead of ``n``.
"""
function localmoran(x::AbstractVector{T} where T<:Real, W::SpatialWeights; permutations::Int = 9999,
    corrected::Bool = true, rng::AbstractRNG = default_rng())::LocalMoran

    n = length(x)
    z = x .- mean(x)    
    Wz = slag(W, z)
    
    m2::Float64 = sum(z.^2) 
    if corrected
        m2 = m2 ./ (n - 1)
    else
        m2 = m2 ./ n
    end

    # Auxiliar function to calculate Local Moran
    function localmoran_calc(zi::Number, wi::AbstractVector, zneighi::AbstractVector)::Float64
        return (zi / m2) .* sum(wi .* zneighi)
    end
    
    # Local Moran
    # I = (z / m2) .* Wz
    I = zeros(n)
    for i in 1:n
        I[i] = localmoran_calc(z[i], weights(W, i), z[neighbors(W, i)]) 
    end

    # Conditional randomizatoin
    Iperms = crand_local(permutations, z, W, localmoran_calc, rng)
    
    larger = sum(Iperms .>= repeat(I, 1, permutations), dims = 2)
    low = (permutations .- larger) .< larger
    larger[low] .= permutations .- larger[low]
    p = (larger .+ 1) ./ (permutations + 1)
    p = vec(p)

    Ipermsstd = vec(std(Iperms, dims = 2, corrected = false))
    Ipermsmean = vec(mean(Iperms, dims = 2))
    zval = (I .- Ipermsmean)  ./ Ipermsstd
    zval = vec(zval)

    # Classification
    q = Array{Symbol}(undef, n)
    for i in 1:n
        if z[i] > 0
            if Wz[i] > 0
                q[i] = :HH
            else
                q[i] = :HL
            end
        else
            if Wz[i] > 0
                q[i] = :LH
            else
                q[i] = :LL
            end
        end
    end

    return LocalMoran(n, I, p, Iperms, Ipermsmean, Ipermsstd, zval, q)

end

score(x::LocalMoran) = x.I;

scoreperms(x::LocalMoran) = x.Iperms;

mean(x::LocalMoran) = x.Ipermsmean;

std(x::LocalMoran) = x.Ipermsstd;

zscore(x::LocalMoran) = x.z;

pvalue(x::LocalMoran) = x.p;

assignments(x::LocalMoran) = x.q;

testname(::LocalMoran) = "Local Moran";

labelsorder(::LocalMoran) = [:ns; :HH; :LL; :LH; :HL];

labelsnames(::LocalMoran) = ["Not Significant"; "High-High"; "Low-Low"; "Low-High"; "High-Low"]

function labelcolor(::LocalMoran, x::String)

    if x == "HH"
        catcolor = :red
        alpha = 1
    elseif x == "LL"
        catcolor = :blue
        alpha = 1
    elseif x == "LH"
        catcolor = :blue
        alpha = 0.4
    elseif x == "HL"
        catcolor = :red
        alpha = 0.4 
    elseif x == "ns"
        catcolor = :lightgrey
        alpha = 0.4
    end

    return (catcolor, alpha)
end