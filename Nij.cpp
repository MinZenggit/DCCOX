// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>

#define EIGEN_NO_DEBUG
#define EIGEN_DONT_PARALLELIZE

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// [[Rcpp::export]]
NumericMatrix NijT(NumericMatrix Trail, double t, double h, int n, int nn, int ord) {
  
  MatrixXd TrailC = MatrixXd::Map(Trail.begin(), nn, 3);
  
  double bw1 = t - 10*h > 0 ? t - 10*h : 0;
  double bw2 = t + 10*h < 1 ? t + 10*h : 1;
  double diff(0.0);
  int i(0), j(0);
  
  MatrixXd NijT = MatrixXd::Zero(n, n);  
  
  static const double inv_sqrt_2pi = 0.3989422804014327;
  
  
  for (int z = 0; z < nn; z++) {
    if (TrailC(z, 2) > bw1 && TrailC(z, 2) < bw2) {
      i = TrailC(z, 0);
      j = TrailC(z, 1);
      diff = (TrailC(z, 2) - t) / h;
      NijT(i - 1, j - 1) += std::pow(inv_sqrt_2pi / h * std::exp(-0.5 * diff * diff), ord);
    }
  }
  
  return Rcpp::wrap(NijT);
}
  