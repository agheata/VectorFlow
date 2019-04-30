# vectorflow
C++ adapter API for integrating vectorized components in a scalar workflow.

Many FLOP-intensive algorithms may profit from the vector pipelines of modern processors; they don’t because they don’t have vectorizable inner loops. The project idea is to implement and benchmark a generic vector flow service that can non-intrusively integrate with arbitrary data processing frameworks and can expose algorithms to the higher-level event loop of these frameworks.

* The service allows to filter from the main data stream to extract data of interest
* Scalar data of arbitrary type can be accumulated in in vectors of a given size
* Once the vectors are filled to the requested size, the scalar data can be transformed in form of a structure of arrays
* The transformed data can be directly fed into a vector-aware implementation of the given algorithm. There is no constraint by the framework on how to implement vectorization: the . Vectorization will happen on elements of the SOA. Vectorization can be compiler driven, but the framework integrates with [VecCore](https://github.com/root-project/veccore) services.
* The output data can be scattered in scalar form and re-integrated in the framework data flow. 

Depending on the intrinsic algorithm gain from SIMD vectorization and better data caching, the overheads introduced by the extra data transformations can be much smaller than the benefits.

The library aims to provide an API and a set of examples demonstrating the transformation of the scalar algorithm in a vectorized one.
