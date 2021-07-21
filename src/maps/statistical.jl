# Statistical classification 
struct BoxPlot <: AbstractStatisticalMapClassificator 
    k::Int
    h::Float64
    BoxPlot(h::Float64 = 1.5) = return new(6, h);
end

struct StdMean <: AbstractStatisticalMapClassificator
    k::Int
    StdMean() = return new(6);
end

struct Percentiles <: AbstractStatisticalMapClassificator
    k::Int
    p::Union{Vector{Int}, Vector{Float64}}
    Percentiles(p::Union{Vector{Int}, Vector{Float64}} = [1, 10, 50, 90, 99]) = return new(length(p) + 1, p);
end

# Bounds
function mapclassifybounds(mcr::BoxPlot, x::Vector{T} where T<:Real)::Tuple{Vector{Float64}, Vector{Float64}}
    h = mcr.h
    q = quantile(x, [0.25, 0.5, 0.75])
    q25 = q[1]
    q50 = q[2]
    q75 = q[3]
    IQR = q75 - q25

    lbound = [-Inf, 
              q25 - h * IQR,
              q25,
              q50,
              q75,
              q75 + h * IQR]
    ubound = [q25 - h * IQR,
              q25,
              q50,
              q75,
              q75 + h * IQR,
              Inf]
    
    lbound, ubound
end

function mapclassifybounds(::StdMean, x::Vector{T} where T<:Real)::Tuple{Vector{Float64}, Vector{Float64}}
    m = mean(x)
    sd = std(x, corrected = true)

    lbounds = [-Inf, -2, -1, 0, 1, 2] * sd .+ m
    ubounds = [-2, -1, 0, 1, 2, Inf]  * sd .+ m

    lbounds, ubounds
end

function mapclassifybounds(mcr::Percentiles, x::Vector{T} where T<:Real)::Tuple{Vector{Float64}, Vector{Float64}}
    p = percentile(x, mcr.p)

    lbound = vcat(minimum(x), p) 
    ubound = vcat(p, maximum(x))

    lbound, ubound
end
