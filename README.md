# Fortran Routines for Least Squares Problems

This package contains the Fortran 77 and Fortran 90 codes accompanying the SIAM Publications printing of "Solving Least Squares Problems" by C. Lawson and R. Hanson [1]. The original routines (most dating back to 1974) available from [netlib](https://www.netlib.org/lawson-hanson/all) are:

* `bndacc` and `bndsol` implement the band-limited sequential algorithm
* `hfti` solves a least squares problem by Householder transformations
* `ldp` solves the least distance programming problem
* `nnls` solves a least squares problem, subject to all variables being nonnegative
* `qrbd` computes the singular value decomposition of a bidiagonal matrix
* `bvls` solves a least squares problem, subject to all variables having upper and lower bounds
* `sva` implements singular value analysis and Levenberg-Marquardt analysis
* `svdrs` computes the singular value decomposition

Some additional utility routines for performing Householder transformations (`h12`), orthogonal rotations (`g1` and `g2`), and generating random integer sequences (`gen`) are also available. More details concerning the routines can be found in the original work [[1]](https://doi.org/10.1137/1.9781611971217).

## Modern Fortran interface

Due to the limitations of early Fortran dialects, the original routines are unwieldy to use, requiring the user to provide several dimensioning variables and manually allocate scratch space. Nevertheless, the usefulness of these routines has prompted users to rewrite them in [C](https://github.com/mutantturkey/nnls_solver), [C++](https://github.com/hmatuschek/eigen3-nnls), [Python](https://github.com/stefanopalmieri/lsqnonneg/blob/master/lsqnonneg.py) and recently [Julia](https://github.com/rdeits/NNLS.jl); or, develop wrappers for languages like [R](https://cran.r-project.org/web/packages/nnls/index.html), [Ruby (through f2c)](https://github.com/mlapshin/nnls) and [Python (SciPy)](https://docs.scipy.org/doc/scipy/reference/generated/scipy.optimize.nnls.html#scipy.optimize.nnls). Ironically, only Fortran developers are stuck with the old legacy interface.

To improve the state of things, this package aims to provide updated versions of the original codes. In the meanwhile, we have prepared some simplified interfaces:

```Fortran
call nnls(A,b,x[,rnorm,ierr])

call bvls(A,b,bnd,x[,rnorm,ierr])

call ldp(G,h,x[,xnorm,ierr])

call hfti(A,b,x,tau[,krank,rnorm])
```

## References

[[1] Lawson, Charles L., and Richard J. Hanson. *Solving least squares problems*. Vol. 15. Siam, 1995.](https://doi.org/10.1137/1.9781611971217)
