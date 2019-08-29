# VECTORFLOW
C++ open source adapter that integrates SIMD vectorized components in a scalar workflow.

Many FLOP-intensive high-energy physics algorithms may profit from the vector pipelines of modern processors, they don’t because they don’t have vectorizable inner loops. The API provides a generic *vector flow* service that can non-intrusively integrate with arbitrary data processing frameworks.

* The service allows the user to work only on data of interest. This is, a *FILTER* step to avoid unnecessary computations.
* The scalar data of arbitrary type can be accumulated in vectors of a given size (defined by the computer architecture). This is, a *GATHER* step to prepare the data for SIMD processing.
* The transformed data, then, can be directly fed into a vector-aware implementation of the given algorithm. This is, a *VECTORIZATION* step, which integrates with [VecCore](https://github.com/root-project/veccore) services.
* At the end, the output data can be pushed back to the initial scalar form. This is, a *SCATTER* step to re-integrate the output into the framework data flow.

Depending on the intrinsic algorithm gain from the SIMD vectorization and better instruction caching, the inherent overheads introduced by the extra data transformations can be much smaller than the benefits. The library aims to provide an API and a set of benchmarked examples demonstrating the transformation of the scalar algorithm into a vectorized one.
