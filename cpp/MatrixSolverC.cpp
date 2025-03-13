// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>

#define EIGEN_NO_DEBUG
#define EIGEN_DONT_PARALLELIZE

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// [[Rcpp::export]]
NumericVector MatrixSolverC(NumericMatrix Ar, NumericVector Br) {
  int n = Ar.nrow();
  
  VectorXd X;
  MatrixXd A = MatrixXd::Map(Ar.begin(), n, n);
  MatrixXd B = VectorXd::Map(Br.begin(), n);
  
  X = A.partialPivLu().solve(B);
  return Rcpp::wrap(X);
}