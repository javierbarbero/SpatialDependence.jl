# Test spatial weights matrix
@testset "Choropleth Maps" begin

    guerry = sdataset("Guerry")
    boston = sdataset("Bostonhsg")

    @testset "EqualIntervals" begin
        mc = mapclassify(EqualIntervals(5), guerry.Litercy)

        @test length(mc) == 5
        @test counts(mc) == [20, 20, 21, 12, 12]
        lbound, ubound = bounds(mc)
        @test lbound == [12, 24.4, 36.8, 49.2, 61.6]
        @test ubound == [24.4, 36.8, 49.2, 61.6, 74]
        @test mc.lower == :open
        @test mc.upper == :closed
        @test levels(mc) == string.(1:5)
        @test assignments(mc)[1:3] == [3, 4, 1]
        @test assignments(mc)[83:85] == [1, 5, 3]
    end

    @testset "Quantiles" begin
        mc = mapclassify(Quantiles(5), guerry.Litercy)

        @test length(mc) == 5
        @test counts(mc) == [18, 19, 15, 16, 17]
        lbound, ubound = bounds(mc)
        @test lbound ≈ [12, 23, 31, 43, 54.4]
        @test ubound ≈ [23, 31, 43, 54.4, 74]
        @test mc.lower == :open
        @test mc.upper == :closed
    end

    @testset "CustomBreaks" begin
        mc = mapclassify(CustomBreaks([20, 30, 60]), guerry.Litercy)

        @test length(mc) == 4
        @test counts(mc) == [14, 18, 41, 12]
        lbound, ubound = bounds(mc)
        @test lbound ≈ [12, 20, 30, 60]
        @test ubound ≈ [20, 30, 60, 74]
        @test mc.lower == :open
        @test mc.upper == :closed
    end

    @testset "NaturalBreaks" begin
        mc = mapclassify(NaturalBreaks(5), guerry.Litercy)

        @test length(mc) == 5
        @test counts(mc) == [15, 24, 22, 15, 9]
        lbound, ubound = bounds(mc)
        @test lbound == [12, 21, 34, 49, 63.0]
        @test ubound == [21, 34, 49, 63, 74.0]
        @test mc.lower == :open
        @test mc.upper == :closed
    end

    @testset "BoxPlot" begin
        mc = mapclassify(BoxPlot(), guerry.Litercy)

        @test length(mc) == 6
        @test counts(mc) == [0, 23, 21, 20 ,21, 0]
        lbound, ubound = bounds(mc)
        @test lbound == [-Inf, -15.5, 25, 38, 52, 92.5]
        @test ubound == [-15.5, 25, 38, 52, 92.5, Inf]
        @test mc.lower == :open
        @test mc.upper == :closed
    end

    @testset "StdMean" begin
        mc = mapclassify(StdMean(), guerry.Litercy)

        @test length(mc) == 6
        @test counts(mc) == [0, 15, 30, 24, 16, 0]
        lbound, ubound = bounds(mc)
        @test lbound ≈ [-Inf, 4.272507, 21.706842, 39.141176, 56.575511, 74.009846] atol = 1e-5
        @test ubound ≈ [4.272507, 21.706842, 39.141176, 56.575511, 74.009846, Inf] atol = 1e-5
        @test mc.lower == :open
        @test mc.upper == :closed
    end

    @testset "Percentiles" begin
        mc = mapclassify(Percentiles(), guerry.Litercy)

        @test length(mc) == 6
        @test counts(mc) == [1, 9, 34, 32, 8, 1]
        lbound, ubound = bounds(mc)
        @test lbound ≈ [12, 12.84, 18, 38, 65.4, 73.16] atol = 1e-5
        @test ubound ≈ [12.84, 18, 38, 65.4, 73.16, 74] atol = 1e-5
        @test mc.lower == :open
        @test mc.upper == :closed
    end

    @testset "Closed-Open" begin
        mc = mapclassify(Quantiles(5), guerry.Litercy, lower = :closed, upper = :open)

        @test length(mc) == 5
        @test counts(mc) == [15, 17, 18, 18, 17]
        lbound, ubound = bounds(mc)
        @test lbound ≈ [12, 23, 31, 43, 54.4]
        @test ubound ≈ [23, 31, 43, 54.4, 74]
        @test mc.lower == :closed
        @test mc.upper == :open
    end

    @testset "Graduated Labels" begin
        mc = mapclassify(Quantiles(5), guerry.Litercy)
        @test maplabels(mc) == 
            ["[12.0, 23.0] (18)"
            "(23.0, 31.0] (19)"
            "(31.0, 43.0] (15)"
            "(43.0, 54.4] (16)"
            "(54.4, 74.0] (17)"]

        @test maplabels(mc, sep = " -- ") == 
            ["[12.0 -- 23.0] (18)"
            "(23.0 -- 31.0] (19)"
            "(31.0 -- 43.0] (15)"
            "(43.0 -- 54.4] (16)"
            "(54.4 -- 74.0] (17)"]

        @test maplabels(mc, counts = false) == 
            ["[12.0, 23.0]"
            "(23.0, 31.0]"
            "(31.0, 43.0]"
            "(43.0, 54.4]"
            "(54.4, 74.0]"]

        mc = mapclassify(StdMean(), guerry.Litercy)
        @test maplabels(mc, digits = 4) == 
            ["[-Inf, 4.2725] (0)"
            "(4.2725, 21.7068] (15)"
            "(21.7068, 39.1412] (30)"
            "(39.1412, 56.5755] (24)"
            "(56.5755, 74.0098] (16)"
            "(74.0098, Inf] (0)"]

        mc = mapclassify(Quantiles(5), guerry.Litercy, lower = :closed, upper = :open)
        @test maplabels(mc) == 
            ["[12.0, 23.0) (15)"
            "[23.0, 31.0) (17)"
            "[31.0, 43.0) (18)"
            "[43.0, 54.4) (18)"
            "[54.4, 74.0] (17)"]
    end

    @testset "Unique" begin
        mc = mapclassify(Unique(), guerry.Region)

        @test length(mc) == 5
        @test counts(mc) == [17, 17, 17, 17, 17]
        @test levels(mc) == ["E", "N", "C", "S", "W"]
    end

    @testset "Unique Labels" begin
        mc = mapclassify(Unique(), guerry.Region)
        @test maplabels(mc) == 
            ["E (17)"
            "N (17)"
            "C (17)"
            "S (17)"
            "W (17)"]

        @test maplabels(mc, counts = false) == 
            ["E"
            "N"
            "C"
            "S"
            "W"]
    end

    @testset "CoLocation" begin
        mc = mapclassify(CoLocation(StdMean()), [guerry.Litercy, guerry.Donatns])

        @test length(mc) == 7
        @test counts(mc) == [61, 0, 0, 19, 3, 2, 0]
        @test levels(mc) == ["Other", "1", "2", "3", "4", "5", "6"]

        @test maplabels(mc) == 
            ["Other (61)"
            "1 (0)"
            "2 (0)"
            "3 (19)"
            "4 (3)"
            "5 (2)"
            "6 (0)"]
        @test maplabels(mc, counts = false) ==
            ["Other"
            "1"
            "2"
            "3"
            "4"
            "5"
            "6"]
    end

    @testset "CoLocation String" begin
        LitHigh = repeat(["Low"], length(guerry.Litercy))
        LitHigh[guerry.Litercy .> mean(guerry.Litercy)] .= "High"
        DonHigh = repeat(["Low"], length(guerry.Donatns))
        DonHigh[guerry.Donatns .> mean(guerry.Donatns)] .= "High"

        mc = mapclassify(CoLocation(Unique()), [LitHigh, DonHigh])

        @test length(mc) == 3
        @test counts(mc) == [49, 10, 26]
        @test levels(mc) == ["Other", "High", "Low"]

        @test maplabels(mc) == 
            ["Other (49)"
            "High (10)"
            "Low (26)"]
        @test maplabels(mc, counts = false) ==
            ["Other"
            "High"
            "Low"]
    end

    @testset "Plot recipes" begin
        # Polygons
        @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}(), guerry.geometry, guerry.Litercy, EqualIntervals())
        # Points
        @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}(), boston.geometry, boston.MEDV, EqualIntervals())
        # Empty group in legend
        @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}(), boston.geometry, boston.MEDV, StdMean())
        # DataFrame and Symbol
        @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}(), guerry, :Litercy, EqualIntervals())
        # Paired palette with k <= 12
        @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}(), guerry, :Region, Unique())
        # Co-location
        @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}(), guerry.geometry, [guerry.Litercy, guerry.Donatns], CoLocation(StdMean()))
        # DataRrame Co-location
        @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}(), guerry, [:Litercy, :Donatns], CoLocation(StdMean()))
    end

end
