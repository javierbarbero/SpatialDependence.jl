# Test spatial weights matrix
@testset "Spatial Weights" begin

    # Create 3x3 polygons
    polygons = SpatialDependence.LatticePolygon[]
    for i in 1:3
        for j in 1:3
            push!(polygons, SpatialDependence.LatticePolygon([[0.0 + i, 0.0 + j], [1.0 + i, 0.0 + j], [1.0 + i, 1.0 + j], [0.0 + i, 1.0 + j]]) )
        end
    end

    # Rook contiguity
    @testset "Rook Contiguity" begin
        W = polyneigh(polygons, criterion = :Rook)
        @test isa(W, SpatialWeights)

        @test nobs(W) == 9
        @test mean(W) ≈ 2.6666666666666665 atol = 1e-5
        @test minimum(W) == 2
        @test maximum(W) == 4
        @test median(W) == 3.0
        @test nislands(W) == 0
        @test islands(W) == []

        @test cardinalities(W) == [2, 3, 2, 3, 4, 3, 2, 3, 2]

        @test neighbors(W, 1) == [2, 4]
        @test neighbors(W, 2) == [1, 3, 5]
        @test neighbors(W, 3) == [2, 6]
        @test neighbors(W, 4) == [1, 5, 7]
        @test neighbors(W, 5) == [2, 4, 6, 8]
        @test neighbors(W, 6) == [3, 5, 9]
        @test neighbors(W, 7) == [4, 8]
        @test neighbors(W, 8) == [5, 7, 9]
        @test neighbors(W, 9) == [6, 8]

        # Default is row-standardized weights
        @test wtransformation(W) == :row

        @test weights(W, 1) == [0.5, 0.5]
        @test weights(W, 2) == [1/3, 1/3, 1/3]
        @test weights(W, 3) == [0.5, 0.5]
        @test weights(W, 4) == [1/3, 1/3, 1/3]
        @test weights(W, 5) == [1/4, 1/4, 1/4, 1/4]
        @test weights(W, 6) == [1/3, 1/3, 1/3]
        @test weights(W, 7) == [0.5, 0.5]
        @test weights(W, 8) == [1/3, 1/3, 1/3]
        @test weights(W, 9) == [0.5, 0.5]

        # Transform to Binary weights
        wtransform!(W, :binary)

        @test weights(W, 1) == [1, 1]
        @test weights(W, 2) == [1, 1, 1]
        @test weights(W, 3) == [1, 1]
        @test weights(W, 4) == [1, 1, 1]
        @test weights(W, 5) == [1, 1, 1, 1]
        @test weights(W, 6) == [1, 1, 1]
        @test weights(W, 7) == [1, 1]
        @test weights(W, 8) == [1, 1, 1]
        @test weights(W, 9) == [1, 1]

        # Copy transformed to Row-standardization
        Wrow = wtransform(W, :row)

        @test wtransformation(Wrow) == :row

        @test weights(Wrow, 1) == [0.5, 0.5]
        @test weights(Wrow, 2) == [1/3, 1/3, 1/3]
        @test weights(Wrow, 3) == [0.5, 0.5]
        @test weights(Wrow, 4) == [1/3, 1/3, 1/3]
        @test weights(Wrow, 5) == [1/4, 1/4, 1/4, 1/4]
        @test weights(Wrow, 6) == [1/3, 1/3, 1/3]
        @test weights(Wrow, 7) == [0.5, 0.5]
        @test weights(Wrow, 8) == [1/3, 1/3, 1/3]
        @test weights(Wrow, 9) == [0.5, 0.5]

        # In-place row standardization
        wtransform!(W, :row)

        @test wtransformation(W) == :row

        @test weights(W, 1) == [0.5, 0.5]
        @test weights(W, 2) == [1/3, 1/3, 1/3]
        @test weights(W, 3) == [0.5, 0.5]
        @test weights(W, 4) == [1/3, 1/3, 1/3]
        @test weights(W, 5) == [1/4, 1/4, 1/4, 1/4]
        @test weights(W, 6) == [1/3, 1/3, 1/3]
        @test weights(W, 7) == [0.5, 0.5]
        @test weights(W, 8) == [1/3, 1/3, 1/3]
        @test weights(W, 9) == [0.5, 0.5]
    end

    # Queen contiguity
    @testset "Queen Contiguity" begin
        W = polyneigh(polygons, criterion = :Queen)
        @test isa(W, SpatialWeights)

        @test nobs(W) == 9
        @test mean(W) ≈ 4.444444444444445 atol = 1e-5
        @test minimum(W) == 3
        @test maximum(W) == 8
        @test median(W) == 5.0
        @test nislands(W) == 0
        @test islands(W) == []

        @test cardinalities(W) == [3, 5, 3, 5, 8, 5, 3, 5, 3]

        @test neighbors(W, 1) == [2, 4, 5]
        @test neighbors(W, 2) == [1, 3, 4, 5, 6]
        @test neighbors(W, 3) == [2, 5, 6]
        @test neighbors(W, 4) == [1, 2, 5, 7, 8]
        @test neighbors(W, 5) == [1, 2, 3, 4, 6, 7, 8, 9]
        @test neighbors(W, 6) == [2, 3, 5, 8, 9]
        @test neighbors(W, 7) == [4, 5, 8]
        @test neighbors(W, 8) == [4, 5, 6, 7, 9]
        @test neighbors(W, 9) == [5, 6, 8]

        # Row-standardized weights
        @test wtransformation(W) == :row

        @test weights(W, 1) == [1/3, 1/3, 1/3]
        @test weights(W, 2) == [1/5, 1/5, 1/5, 1/5, 1/5]
        @test weights(W, 3) == [1/3, 1/3, 1/3]
        @test weights(W, 4) == [1/5, 1/5, 1/5, 1/5, 1/5]
        @test weights(W, 5) == [1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8, 1/8]
        @test weights(W, 6) == [1/5, 1/5, 1/5, 1/5, 1/5]
        @test weights(W, 7) == [1/3, 1/3, 1/3]
        @test weights(W, 8) == [1/5, 1/5, 1/5, 1/5, 1/5]
        @test weights(W, 9) == [1/3, 1/3, 1/3]        

        # Binary weights
        wtransform!(W, :binary)
        @test wtransformation(W) == :binary

        @test weights(W, 1) == [1, 1, 1]
        @test weights(W, 2) == [1, 1, 1, 1, 1]
        @test weights(W, 3) == [1, 1, 1]
        @test weights(W, 4) == [1, 1, 1, 1, 1]
        @test weights(W, 5) == [1, 1, 1, 1, 1, 1, 1, 1]
        @test weights(W, 6) == [1, 1, 1, 1, 1]
        @test weights(W, 7) == [1, 1, 1]
        @test weights(W, 8) == [1, 1, 1, 1, 1]
        @test weights(W, 9) == [1, 1, 1]

        @test_nowarn show(IOBuffer(), W)
    end

    # Test Multipolygon with multipolygon 1 at the left and at the right
    @testset "Multipolygon" begin
        mpolygons = SpatialDependence.LatticeMultiPolygon[]
        push!(mpolygons,       
            SpatialDependence.LatticeMultiPolygon([
                    SpatialDependence.LatticePolygon([[1.0, 0.0], [2.0, 0.0], [2.0, 1.0], [1.0, 1.0]]),
                    SpatialDependence.LatticePolygon([[5.0, 0.0], [6.0, 0.0], [6.0, 1.0], [5.0, 1.0]])])
            )

        for i in 2:4
            push!(mpolygons, SpatialDependence.LatticeMultiPolygon(SpatialDependence.LatticePolygon([[0.0 + i, 0.0], [1.0 + i, 0.0 ], [1.0 + i, 1.0], [0.0 + i, 1.0]])) )
        end

        W = polyneigh(mpolygons) 

        @test neighbors(W, 1) == [2, 4]
        @test neighbors(W, 2) == [1, 3]
        @test neighbors(W, 3) == [2, 4]
        @test neighbors(W, 4) == [1, 3]
    end

    # Test tolerance for edge detection
    @testset "Polygon tolerance" begin
        polygonstol = SpatialDependence.SpatialDependence.LatticePolygon[]
        push!(polygonstol, SpatialDependence.LatticePolygon([[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0]]) )
        push!(polygonstol, SpatialDependence.LatticePolygon([[1.0 + 0.25, 0.0], [2.0 + 0.25, 0.0], [2.0 + 0.25, 1.0], [1.0 + 0.25, 1.0]]) )
        push!(polygonstol, SpatialDependence.LatticePolygon([[2.0 + 0.50, 0.0], [3.0 + 0.50, 0.0], [3.0 + 0.50, 1.0], [2.0 + 0.50, 1.0]]) )

        W = polyneigh(polygonstol) 

        @test neighbors(W, 1) == []
        @test neighbors(W, 2) == []
        @test neighbors(W, 3) == []

        @test nislands(W) == 3
        @test islands(W) == [1, 2, 3]

        W = polyneigh(polygonstol, tol = 0.25) 

        @test neighbors(W, 1) == [2]
        @test neighbors(W, 2) == [1, 3]
        @test neighbors(W, 3) == [2]

        @test nislands(W) == 0
        @test islands(W) == []
    end

    # Generate points for test
    points = SpatialDependence.LatticePoint[]
    for i in 1:3
        for j in 1:3
            push!(points, SpatialDependence.LatticePoint([0.0 + i, 0.0 + j]))
        end
    end
    
    # Points distance threshold
    @testset "Distance Threshold" begin
        W = dnearneigh(points, threshold = 1.5)
        @test isa(W, SpatialWeights)

        @test nobs(W) == 9
        @test mean(W) ≈ 4.444444444444445 atol = 1e-5
        @test minimum(W) == 3
        @test maximum(W) == 8
        @test median(W) == 5.0
        @test nislands(W) == 0
        @test islands(W) == []

        @test neighbors(W, 1) == [2, 4, 5]
        @test neighbors(W, 2) == [1, 3, 4, 5, 6]
        @test neighbors(W, 3) == [2, 5, 6]
        @test neighbors(W, 4) == [1, 2, 5, 7, 8]
        @test neighbors(W, 5) == [1, 2, 3, 4, 6, 7, 8, 9]
        @test neighbors(W, 6) == [2, 3, 5, 8, 9]
        @test neighbors(W, 7) == [4, 5, 8]
        @test neighbors(W, 8) == [4, 5, 6, 7, 9]
        @test neighbors(W, 9) == [5, 6, 8]
    end

    @testset "Polygon Missing" begin
        PM = Vector{Union{Missing,SpatialDependence.LatticePolygon}}(undef, 10)
        PM[1:9] = polygons[1:9]
        PM[10] = missing

        @test typeof(PM) == Vector{Union{Missing,SpatialDependence.LatticePolygon}}
        @test_throws ArgumentError polyneigh(PM) # Missin geometry error
    end

    # Points k near neighbors
    @testset "K near neighbors " begin
        W = knearneigh(points, k = 3)
        @test isa(W, SpatialWeights)

        @test nobs(W) == 9
        @test mean(W) ≈ 3 atol = 1e-5
        @test minimum(W) == 3
        @test maximum(W) == 3
        @test median(W) == 3.0
        @test nislands(W) == 0
        @test islands(W) == []

        @test neighbors(W, 1) == [2, 4, 5]
        @test neighbors(W, 2) == [1, 3, 5]
        @test neighbors(W, 3) == [2, 5, 6]
        @test neighbors(W, 4) == [1, 5, 7]
        @test neighbors(W, 5) == [2, 4, 8]
        @test neighbors(W, 6) == [3, 5, 9]
        @test neighbors(W, 7) == [4, 5, 8]
        @test neighbors(W, 8) == [5, 7, 9]
        @test neighbors(W, 9) == [5, 6, 8]
    end

    # Polygon mean center
    @testset "Polygon meancenter" begin
        x, y = meancenter(reggeomlattice(3, 3, direction = :rightup))
        @test x == [1.5, 2.5, 3.5, 1.5, 2.5, 3.5, 1.5, 2.5, 3.5]
        @test y == [1.5, 1.5, 1.5, 2.5, 2.5, 2.5, 3.5, 3.5, 3.5]
    end

    # Polygon centroids
    @testset "Polygon centroid" begin
        x, y = centroid(reggeomlattice(3, 3, direction = :rightup))
        @test x == [1.5, 2.5, 3.5, 1.5, 2.5, 3.5, 1.5, 2.5, 3.5]
        @test y == [1.5, 1.5, 1.5, 2.5, 2.5, 2.5, 3.5, 3.5, 3.5]
    end

    # Spatial Weights from Matrix
    @testset "From Matrix " begin
        W = SpatialWeights([0 1 0; 1 0 1; 0 1 0])
        @test isa(W, SpatialWeights)

        @test nobs(W) == 3
        @test mean(W) ≈ 4/3 atol = 1e-5
        @test minimum(W) == 1
        @test maximum(W) == 2
        @test median(W) == 1.0
        @test nislands(W) == 0
        @test islands(W) == []

        @test neighbors(W, 1) == [2]
        @test neighbors(W, 2) == [1, 3]
        @test neighbors(W, 3) == [2]

        # Row-standardize by default
        @test wtransformation(W) == :row
        @test weights(W, 1) == [1]
        @test weights(W, 2) == [0.5, 0.5]
        @test weights(W, 3) == [1]

        # Do not standardize
        W = SpatialWeights([0 2 0; 3 0 4; 0 5 0], standardize = false)
        @test wtransformation(W) == :original
        @test weights(W, 1) == [2]
        @test weights(W, 2) == [3, 4]
        @test weights(W, 3) == [5] 

    end

    # Test conversions
    @testset "Matrix converion" begin   
        W = polyneigh(polygons, criterion = :Rook)

        @test Matrix(W) ==
        [
            0.0   0.5   0.0   0.5   0.0   0.0   0.0   0.0   0.0
            1/3   0.0   1/3   0.0   1/3   0.0   0.0   0.0   0.0
            0.0   0.5   0.0   0.0   0.0   0.5   0.0   0.0   0.0
            1/3   0.0   0.0   0.0   1/3   0.0   1/3   0.0   0.0
            0.0   0.25  0.0   0.25  0.0   0.25  0.0   0.25  0.0
            0.0   0.0   1/3   0.0   1/3   0.0   0.0   0.0   1/3
            0.0   0.0   0.0   0.5   0.0   0.0   0.0   0.5   0.0
            0.0   0.0   0.0   0.0   1/3   0.0   1/3   0.0   1/3
            0.0   0.0   0.0   0.0   0.0   0.5   0.0   0.5   0.0
        ]

        # Binary weights
        wtransform!(W, :binary)
        @test Matrix(W) == 
        [
            0.0  1.0  0.0  1.0  0.0  0.0  0.0  0.0  0.0;
            1.0  0.0  1.0  0.0  1.0  0.0  0.0  0.0  0.0;
            0.0  1.0  0.0  0.0  0.0  1.0  0.0  0.0  0.0;
            1.0  0.0  0.0  0.0  1.0  0.0  1.0  0.0  0.0;
            0.0  1.0  0.0  1.0  0.0  1.0  0.0  1.0  0.0;
            0.0  0.0  1.0  0.0  1.0  0.0  0.0  0.0  1.0;
            0.0  0.0  0.0  1.0  0.0  0.0  0.0  1.0  0.0;
            0.0  0.0  0.0  0.0  1.0  0.0  1.0  0.0  1.0;
            0.0  0.0  0.0  0.0  0.0  1.0  0.0  1.0  0.0 
        ]

        @test Matrix(sparse(W)) == Matrix(W)        

    end

    # Regular lattice
    @testset "Regular lattice" begin   
        regrightup = reggeomlattice(3, 3, direction = :rightup)
        regrightdown = reggeomlattice(3, 3, direction = :rightdown)
        regleftdown = reggeomlattice(3, 3, direction = :leftdown)
        regleftup = reggeomlattice(3, 3, direction = :leftup)

        @test GeoInterface.coordinates(regrightup[1])   == [[[[1.0, 1.0], [2.0, 1.0], [2.0, 2.0], [1.0, 2.0]]]]
        @test GeoInterface.coordinates(regrightdown[1]) == [[[[1.0, 3.0], [2.0, 3.0], [2.0, 4.0], [1.0, 4.0]]]]
        @test GeoInterface.coordinates(regleftdown[1])  == [[[[3.0, 3.0], [4.0, 3.0], [4.0, 4.0], [3.0, 4.0]]]]
        @test GeoInterface.coordinates(regleftup[1])    == [[[[3.0, 1.0], [4.0, 1.0], [4.0, 2.0], [3.0, 2.0]]]]

        WRook = polyneigh(regrightup, criterion = :Rook)
        @test minimum(WRook) == 2
        @test maximum(WRook) == 4
        @test median(WRook) == 3.0

        WhtQueen = polyneigh(regrightup, criterion = :Queen)
        @test minimum(WhtQueen) == 3
        @test maximum(WhtQueen) == 8
        @test median(WhtQueen) == 5.0

        regpoints = reggeomlattice(3, 3, gtype = :point)
        Wpoint = dnearneigh(regpoints, threshold = 1.0)
        @test minimum(Wpoint) == 2
        @test maximum(Wpoint) == 4
        @test median(WRook) == 3.0
    end

end
