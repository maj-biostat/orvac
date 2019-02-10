

run_simulation <- function(){

  # involves the rcpp parallel routines

  # m <- matrix(as.numeric(c(1:10000000)), nrow = 1000, ncol = 1000)
  #
  # # ensure that serial and parallel versions give the same result
  # stopifnot(identical(matrixSqrt(m), parallelMatrixSqrt(m)))
  #
  # # compare performance of serial and parallel
  # library(rbenchmark)
  # res <- benchmark(matrixSqrt(m),
  #                  parallelMatrixSqrt(m),
  #                  order="relative")
  # res[,1:4]
}


