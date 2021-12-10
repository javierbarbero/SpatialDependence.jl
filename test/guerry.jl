# Test using the Guerry's Moral statistics of France dataset
@testset "Guerry" begin
    
    guerry = sdataset("Guerry")
    W = polyneigh(guerry)

    # In this dataset there is no difference between Queen or Rook contiguity
    @testset "W contiguity" begin
        @test nobs(W) == 85
        @test mean(W) ≈ 4.941176 atol = 1e-5
        @test median(W) ≈ 5.0
        @test minimum(W) == 2
        @test maximum(W) == 8
        @test nislands(W) == 0
        @test islands(W) == []

        # Connectivity histogram
        @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}(), W)
    end

    @testset "Spatial lag" begin   
        @test W * guerry.Litercy == slag(W, guerry.Litercy)
    end

    @testset "Moran's I" begin
        mguerry = moran(guerry.Litercy, W, rng = StableRNG(1234567), permutations = 9999)

        @test score(mguerry) ≈ 0.717605 atol = 1e-5    
        @test length(scoreperms(mguerry)) == 9999
        @test mean(mguerry) ≈ -0.012594 atol = 1e-5
        @test std(mguerry) ≈ 0.070790 atol = 1e-5
        @test zscore(mguerry) ≈ 10.315064 atol = 1e-5
        @test pvalue(mguerry) ≈ 0.0001 atol = 1e-5   
        @test SpatialDependence.expected(mguerry) == - 1 / (length(guerry.Litercy) - 1)
        @test SpatialDependence.testname(mguerry) == "Moran's I"

        @test_nowarn show(IOBuffer(), mguerry)

        # Test Moran's I distribution
        @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}(), mguerry)
    end

    @testset "Geary's c" begin
        cguerry = geary(guerry.Litercy, W, rng = StableRNG(1234567), permutations = 9999)

        @test score(cguerry) ≈ 0.250202 atol = 1e-5
        @test length(scoreperms(cguerry)) == 9999
        @test mean(cguerry) ≈ 1.000411 atol = 1e-5
        @test std(cguerry) ≈ 0.072280 atol = 1e-5
        @test zscore(cguerry) ≈ -10.379150 atol = 1e-5
        @test pvalue(cguerry) ≈ 0.0001 atol = 1e-5
        @test SpatialDependence.expected(cguerry) == 1
        @test SpatialDependence.testname(cguerry) == "Geary's c"

        @test_nowarn show(IOBuffer(), cguerry)

        # Test Geary's c distribution
        @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}(), cguerry)
    end

    @testset "Polygon mean center" begin
        cx, cy = meancenter(guerry.geometry)
        @test cx[1] ≈ 836732.867469879565760 atol = 1e-5
        @test cy[1] ≈ 2125844.293172690551728 atol = 1e-5
        @test cx[50] ≈ 814647.356215213309042 atol = 1e-5
        @test cy[50] ≈ 2350156.217068645637482 atol = 1e-5
    end

    @testset "Polygon centroid" begin
        cx, cy = centroid(guerry.geometry)
        @test cx[1] ≈ 832852.278780025895685 atol = 1e-5
        @test cy[1] ≈ 2126600.575969800818712 atol = 1e-5
        @test cx[50] ≈ 815193.599669872317463 atol = 1e-5
        @test cy[50] ≈ 2349544.092820507008582 atol = 1e-5
    end

    @testset "Moran 5 KNN" begin
        cx, cy = centroid(guerry.geometry)
        Wk5 = knearneigh(cx, cy, k = 5)
        morank5 = moran(guerry.Litercy, Wk5)
        @test score(morank5) ≈ 0.678857 atol = 1e-5
    end

    @testset "Local Moran" begin
        lmguerry = localmoran(guerry.Litercy, W, rng = StableRNG(1234567), permutations = 9999)

        @test score(lmguerry)[1:5] ≈ [-0.039511; 0.599223; 1.416860; 0.081434; -0.013868] atol = 1e-5
        @test score(lmguerry)[81:85] ≈ [0.270916; 0.603623; 1.488529; 2.145242; 0.203191] atol = 1e-5
        @test size(scoreperms(lmguerry), 1) == 85
        @test size(scoreperms(lmguerry), 2) == 9999
        @test mean(lmguerry)[1:5] ≈ [9.1e-5; -0.006589; -0.025567; -0.002702; -0.029142] atol = 1e-5
        @test mean(lmguerry)[81:85] ≈ [-0.004464; -0.006655; -0.02424; -0.021923; -0.003153] atol = 1e-5
        @test std(lmguerry)[1:5] ≈ [0.060504; 0.270397; 0.588832; 0.193349; 0.96366] atol = 1e-5
        @test std(lmguerry)[81:85] ≈ [0.315696; 0.320665; 0.58715; 0.514386; 0.197625] atol = 1e-5
        @test zscore(lmguerry)[1:5] ≈ [-0.654535; 2.240454; 2.44964; 0.435151; 0.01585] atol = 1e-5
        @test zscore(lmguerry)[81:85] ≈ [0.872291; 1.903161; 2.576463; 4.213106; 1.044118] atol = 1e-5
        @test pvalue(lmguerry)[1:5] ≈ [0.2536; 0.0148; 0.0044; 0.3323; 0.4805] atol = 1e-4   
        @test pvalue(lmguerry)[81:85] ≈ [0.2015; 0.0239; 0.0021; 0.0001; 0.1532] atol = 1e-4   

        @test assignments(lmguerry)[1:5] == [:LH, :HH, :LL, :HH, :HL]
        @test assignments(lmguerry)[81:85] == [:LL, :LL, :LL, :HH, :HH]
        @test count(issignificant(lmguerry, 0.05, adjust = :none)) == 38
        @test count(issignificant(lmguerry, 0.01, adjust = :none)) == 21
        @test count(issignificant(lmguerry, 0.05, adjust = :bonferroni)) == 6
        @test count(issignificant(lmguerry, 0.05, adjust = :fdr)) == 23

        @test SpatialDependence.testname(lmguerry) == "Local Moran"
        @test_nowarn show(IOBuffer(), lmguerry)

        # Test LISA Cluster Map
        @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}(), guerry, lmguerry)
        @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}(), guerry.geometry, lmguerry)

        # Test dividing by n instead of (n - 1)
        lmguerry = localmoran(guerry.Litercy, W, rng = StableRNG(1234567), permutations = 0, corrected = false)
        @test score(lmguerry)[1:5] ≈ [-0.039981; 0.606357; 1.433727; 0.082403; -0.014033] atol = 1e-5
        @test score(lmguerry)[81:85] ≈ [0.274141; 0.610809; 1.506250; 2.170780; 0.205610] atol = 1e-5
    end

end
