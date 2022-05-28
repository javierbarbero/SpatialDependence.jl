# Plot  recipes for choropleth maps

# Map shape coordinates for Plots
function mapshapecoords(gtype::Union{GI.PolygonTrait, GI.MultiPolygonTrait},P)::Tuple{Vector{Float64}, Vector{Float64}}
    scoords = broadcast(GeoInterfaceRecipes._coordvecs, Ref(gtype), P)

    x = map(a -> a[1], scoords)            
    y = map(a -> a[2], scoords)

    x = reduce(vcat, x)
    y = reduce(vcat, y)  

    x, y
end

function mapshapecoords(::GI.PointTrait, P)::Tuple{Vector{Float64}, Vector{Float64}}
    scoords = broadcast(GeoInterfaceRecipes._coordvecs, Ref(GI.PointTrait()), P)

    x = map(a -> a[1][1], scoords)
    y = map(a -> a[1][2], scoords)

    x, y
end

# Choropleth Map
@recipe function f(A::Any, colorvar::Union{AbstractVector,Symbol}, mcr::AbstractMapClassificator)

    if istable(A)
        # If Table or DataFrame
        (:geometry in propertynames(A)) || throw(ArgumentError("table does not have :geometry information"))
        P = A.geometry

        if isa(colorvar, Symbol) 
            cvar = A[!, colorvar]
        elseif isa(colorvar, Vector{Symbol})
            cvar = [A[!, c] for c in colorvar]
        else
            cvar = colorvar
        end
        #(isa(colorvar, Symbol) || isa(colorvar, Vector{Symbol})) ? cvar = A[!, colorvar] : cvar = colorvar
        
    else
        P = A
        cvar = colorvar
    end

    # Check is geometry
    all(GI.isgeometry.(P)) || throw(ArgumentError("Unknown geometry"))

    isCoLocation = false
    if isa(mcr, AbstractCoLocationMapClassificator)
        isCoLocation = true
    else
        length(P) == length(cvar) || throw(ArgumentError("dimensions must match: A has ($(length(P))), colorvar has ($(length(cvar)))"))
    end

    gtype = GI.geomtrait(P[1])

    # Default color palette
    defpalette, defrev = defaultpalette(mcr)

    # Get plot attributes
    lower  = get(plotattributes, :lower, :open)
    upper  = get(plotattributes, :upper, :closed)
    digits = get(plotattributes, :digits, 2)
    sep    = get(plotattributes, :sep, ", ")
    counts = get(plotattributes, :counts, true)
    rev    = get(plotattributes, :rev, defrev)

    # Classify values and get labels
    if isa(mcr, AbstractGraduatedMapClassificator)
        mc = mapclassify(mcr, cvar, lower = lower, upper = upper)
        labels = maplabels(mc, digits = digits, sep = sep, counts = counts)
    elseif isa(mcr, AbstractUniqueMapClassificator)
        mc = mapclassify(mcr, cvar)
        labels = maplabels(mc, counts = counts)
    elseif isa(mcr, AbstractCoLocationMapClassificator)
        mc = mapclassify(mcr, cvar, lower = lower, upper = upper)
        labels = maplabels(mc, counts = counts)
    end

    k = length(mc)
    group = assignments(mc)

    # If palette is :Paired adjust for number of categories
    if defpalette == :Paired
        if k <= 12 defpalette = Symbol(string(defpalette) * "_" * string(k)) end
    end

    # Default attributes for maps
    legend --> true
    grid --> false
    showaxis --> false
    ticks --> false
    aspect_ratio --> :equal
    color_palette --> defpalette

    # Get colors from pallete
    mappal = palette(plotattributes[:color_palette], k - isCoLocation, rev = rev)
    
    # For each category in the map
    for i in 1:k

        # Subset of geometries for the category
        subset = P[group .== i]

        # color for the category
        if isCoLocation
            # Assing grey color for CoLocation when there is no match
            if i > 1
                catcolor = mappal[i - 1]
            else
                catcolor = :lightgrey
            end
        else
            catcolor = mappal[i]
        end
        
        # If empty subset, plot a NaN point to display the text in the legend but not the color
        if length(subset) == 0
            @series begin
                seriestype := :scatter
                seriescolor := catcolor
                label := labels[i]
                [NaN], [NaN]                
            end
            continue
        end

        # Get shape map coordinates and plot polygon or point
        x, y = mapshapecoords(gtype, P[group .== i])

        if isa(gtype, GI.PolygonTrait) || isa(gtype, GI.MultiPolygonTrait)
            @series begin
                seriestype := :shape
                seriescolor := catcolor
                label := labels[i]
                (x, y)
            end
        elseif isa(gtype, GI.PointTrait)
            @series begin
                seriestype := :scatter
                seriescolor := catcolor
                label := labels[i]
                x, y
            end
        end
    end

end

# LISA Cluster Map
@recipe function f(A::Any, lisavar::AbstractLocalSpatialAutocorrelation)

    if istable(A)
        # If Table or DataFrame
        (:geometry in propertynames(A)) || throw(ArgumentError("table does not have :geometry information"))
        P = A.geometry
    else
        P = A
    end

    # Check is geometry
    all(GI.isgeometry.(P)) || throw(ArgumentError("Unknown geometry"))

    gtype = GI.geomtrait(P[1])

    # Get plot attributes
    counts = get(plotattributes, :counts, true)
    pthreshold = get(plotattributes, :sig, 0.05)
    adjust = get(plotattributes, :adjust, :none)

   q = assignments(lisavar, pthreshold, adjust = adjust)

    mc = mapclassify(Unique(labelsorder(lisavar)), q)
    labels = maplabels(mc, counts = counts)

    k = length(mc)
    group = assignments(mc)
    grouplabs = mc.grouplabs

    # Default attributes for maps
    legend --> true
    grid --> false
    showaxis --> false
    ticks --> false
    aspect_ratio --> :equal    

    # For each category in the map
    for i in 1:k

        # Subset of geometries for the category
        subset = P[group .== i]

        if grouplabs[i] == "ns"
            catcolor = :lightgrey
            alpha = 0.4
        end

        catcolor, alpha = labelcolor(lisavar, grouplabs[i] )
        
        # If empty subset, plot a NaN point to display the text in the legend but not the color
        if length(subset) == 0
            @series begin
                seriestype := :scatter
                seriescolor := catcolor
                alpha := alpha
                label := labels[i]
                [NaN], [NaN]                
            end
            continue
        end

        # Get shape map coordinates and plot polygon or point
        x, y = mapshapecoords(gtype, P[group .== i])

        if isa(gtype, GI.PolygonTrait) || isa(gtype, GI.MultiPolygonTrait)
            @series begin
                seriestype := :shape
                seriescolor := catcolor
                alpha := alpha
                label := labels[i]                
                (x, y)
            end
        elseif isa(gtype, GI.PointTrait)
            @series begin
                seriestype := :scatter
                seriescolor := catcolor
                alpha := alpha
                label := labels[i]
                x, y
            end
        end
    end

end
