# Test using the Guerry's Moral statistics of France dataset
@testset "Guerry" begin
    
    guerry = sdataset("Guerry")
    W = polyneigh(guerry)
    Wbin = wtransform(W, :binary)

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

    @testset "Local Geary" begin
        lcguerry = localgeary(guerry.Litercy, W, rng = StableRNG(1234567), permutations = 9999)

        @test score(lcguerry)[1:5] ≈ [1.191784; 0.275807; 0.439757; 0.949973; 3.134225] atol = 1e-5
        @test score(lcguerry)[81:85] ≈ [0.272243; 0.325157; 0.446885; 0.20288; 0.734317] atol = 1e-5
        @test size(scoreperms(lcguerry), 1) == 85
        @test size(scoreperms(lcguerry), 2) == 9999
        @test mean(lcguerry)[1:5] ≈ [1.013359; 1.468334; 3.269937; 1.157147; 3.95656] atol = 1e-5
        @test mean(lcguerry)[81:85] ≈ [1.412199; 1.664964; 3.272061; 2.744478; 1.209385] atol = 1e-5
        @test std(lcguerry)[1:5] ≈ [0.564736; 0.55419; 1.388559; 0.522587; 1.82078] atol = 1e-5
        @test std(lcguerry)[81:85] ≈ [0.950857; 0.883169; 1.389196; 0.968423; 0.490151] atol = 1e-5
        @test zscore(lcguerry)[1:5] ≈ [0.315944; -2.151836; -2.038214; -0.396441; -0.451638] atol = 1e-5
        @test zscore(lcguerry)[81:85] ≈ [-1.198872; -1.517045; -2.033677; -2.624471; -0.969228] atol = 1e-5
        @test pvalue(lcguerry)[1:5] ≈ [0.3506; 0.0036; 0.003; 0.373; 0.3455] atol = 1e-4   
        @test pvalue(lcguerry)[81:85] ≈ [0.0722; 0.0325; 0.0031; 0.0002; 0.1712] atol = 1e-4   

        @test assignments(lcguerry)[1:5] == [:N, :P, :P, :P, :P]
        @test assignments(lcguerry)[81:85] == [:P, :P, :P, :P, :P]
        @test count(issignificant(lcguerry, 0.05, adjust = :none)) == 49
        @test count(issignificant(lcguerry, 0.01, adjust = :none)) == 29
        @test count(issignificant(lcguerry, 0.05, adjust = :bonferroni)) == 9
        @test count(issignificant(lcguerry, 0.05, adjust = :fdr)) == 40

        @test SpatialDependence.testname(lcguerry) == "Local Geary"
        @test_nowarn show(IOBuffer(), lcguerry)

        # Test LISA Cluster Map
        @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}(), guerry, lcguerry)
        @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}(), guerry.geometry, lcguerry)

        # Test :moran categories
        lcguerry = localgeary(guerry.Litercy, W, rng = StableRNG(1234567), permutations = 9999, categories = :moran)
        @test assignments(lcguerry)[1:5] == [:NE, :HH, :LL, :HH, :OP]
        @test assignments(lcguerry)[81:85] == [:LL, :LL, :LL, :HH, :HH]
        @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}(), guerry, lcguerry)

        # Test :moran categories with the Donatn's variables to test the :NE label when plotting maps
        lcguerry = localgeary(guerry.Donatns, W, rng = StableRNG(1234567), permutations = 9999, categories = :moran)
        @test SpatialDependence.labelsnames(lcguerry) == ["Not Significant"; "High-High"; "Low-Low"; "Other Positive"; "Negative"]
        @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}(), guerry, lcguerry)

        # Test dividing by n instead of (n - 1)
        lcguerry = localgeary(guerry.Litercy, W, rng = StableRNG(1234567), permutations = 0, corrected = false)
        @test score(lcguerry)[1:5] ≈ [1.205972; 0.279091; 0.444992; 0.961282; 3.171538] atol = 1e-5
        @test score(lcguerry)[81:85] ≈ [0.275484; 0.329028; 0.452205; 0.205295; 0.743059] atol = 1e-5
    end

    @testset "Getis-Ord" begin
        goguerry = getisord(guerry.Litercy, W, rng = StableRNG(1234567), permutations = 9999, star = false)

        @test score(goguerry)[1:5] ≈ [0.013602; 0.016636; 0.00684; 0.01303; 0.011971] atol = 1e-5
        @test score(goguerry)[81:85] ≈ [0.009624; 0.007924; 0.006588; 0.020725; 0.014329] atol = 1e-5
        @test size(scoreperms(goguerry), 1) == 85
        @test size(scoreperms(goguerry), 2) == 9999
        @test mean(goguerry)[1:5] ≈ [0.011893; 0.011896; 0.011901; 0.011893; 0.011923] atol = 1e-5
        @test mean(goguerry)[81:85] ≈ [0.011901; 0.011897; 0.011896; 0.011899; 0.011896] atol = 1e-5
        @test std(goguerry)[1:5] ≈ [0.002611; 0.002116; 0.002066; 0.002612; 0.003011] atol = 1e-5
        @test std(goguerry)[81:85] ≈ [0.002611; 0.002087; 0.00206; 0.002095; 0.00233] atol = 1e-5
        @test zscore(goguerry)[1:5] ≈ [0.654535; 2.240454; -2.44964; 0.435151; 0.01585] atol = 1e-5
        @test zscore(goguerry)[81:85] ≈ [-0.872291; -1.903161; -2.576463; 4.213106; 1.044118] atol = 1e-5
        @test pvalue(goguerry)[1:5] ≈ [0.2628; 0.0147; 0.0042; 0.3323; 0.4833] atol = 1e-4   
        @test pvalue(goguerry)[81:85] ≈ [0.1931; 0.0224; 0.0018; 0.0001; 0.157] atol = 1e-4   

        @test assignments(goguerry)[1:5] == [:H, :H, :L, :H, :H]
        @test assignments(goguerry)[81:85] == [:L, :L, :L, :H, :H]
        @test count(issignificant(goguerry, 0.05, adjust = :none)) == 38
        @test count(issignificant(goguerry, 0.01, adjust = :none)) == 21
        @test count(issignificant(goguerry, 0.05, adjust = :bonferroni)) == 6
        @test count(issignificant(goguerry, 0.05, adjust = :fdr)) == 23

        @test SpatialDependence.testname(goguerry) == "Getis-Ord Gi"
        @test_nowarn show(IOBuffer(), goguerry)

        # Test LISA Cluster Map
        @test_nowarn RecipesBase.apply_recipe(Dict{Symbol, Any}(), guerry, goguerry)

        # Getis-Ord Gi*
        goguerry = getisord(guerry.Litercy, W, rng = StableRNG(1234567), permutations = 9999, star = true)

        @test score(goguerry)[1:5] ≈ [0.012985; 0.016231; 0.006398; 0.013045; 0.013977] atol = 1e-5
        @test score(goguerry)[81:85] ≈ [0.009318; 0.007815; 0.006183; 0.020095; 0.014127] atol = 1e-5

        @test SpatialDependence.testname(goguerry) == "Getis-Ord Gi*"

        # Getis-Ord Gi* with binary weights
        goguerry = getisord(guerry.Litercy, Wbin, rng = StableRNG(1234567), permutations = 9999, star = true)

        @test score(goguerry)[1:5] ≈ [0.064923; 0.113616; 0.044785; 0.065224; 0.055906] atol = 1e-5
        @test score(goguerry)[81:85] ≈ [0.046589; 0.054704; 0.043282; 0.140667; 0.084761] atol = 1e-5
    end

end
