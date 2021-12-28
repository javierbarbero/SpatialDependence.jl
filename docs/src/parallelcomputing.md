```@meta
CurrentModule = SpatialDependence
```

# Parallel Computing

The package takes advantage of Julia multi-threading features to increase performance if Julia is started with multiple threads.

See official Julia documentation on [starting Julia with multiple threats](https://docs.julialang.org/en/v1/manual/multi-threading/#Starting-Julia-with-multiple-threads): 

The following functions take advantage of multi-threading to increase their performance:

| Function             | Function Name | What is parallelized?                 |
|:---------------------|:--------------|:--------------------------------------|
| Polygon contiguity   | `polyneigh`   | Bouning box overlaps and polygon hits |
| Global Moran's I     | `moran`       | Permutation test                      |
| Global Geary's c     | `geary`       | Permutation test                      |
| Local Moran          | `localmoran`  | Conditional permutation test          |
| Local Geary          | `localgeary`  | Conditional permutation test          |
| Getis-Ord Statistics | `getisord`    | Conditional permutation test          |
