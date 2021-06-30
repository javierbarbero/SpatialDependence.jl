# Recipes for Moran plot

@recipe function f(x::Vector{Float64}, W::SpatialWeights, zstandardize::Bool = false)

    Wx = slag(W, x)

    if zstandardize
        z = standardize(ZScoreTransform, x)
        Wz = standardize(ZScoreTransform, Wx)
    else
        z = x
        Wz = Wx
    end

    title --> "Moran Scatterplot"
    xguide --> "Attribute"
    yguide --> "Spatial Lag"
    legend --> false
    grid --> false

    # Variable and spatial lag of variable
    @series begin
        seriestype := :scatter
        smooth := true
        z, Wz
    end

    # Vertical line at variable mean
    zmean = mean(z)
    @series begin
        seriestype := :vline
        linestyle := :dash
        seriescolor --> :red
        [zmean]
    end

    # Horizontal line at W * variable mean
    Wzmean = mean(Wz)
    @series begin
        seriestype := :hline
        linestyle := :dash
        seriescolor --> :red
        [Wzmean]
    end

end