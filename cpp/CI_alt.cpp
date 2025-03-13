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
Rcpp::List CI(NumericMatrix Trail, NumericVector Zij, NumericVector xe, double t, double h1, double h2, int n, int nn, int p) {
  
  if (Rf_isNull(Zij.attr("dim"))) {
    throw std::runtime_error("'x' does not have 'dim' attibute.");
  }
  Rcpp::Dimension d = Zij.attr("dim");
  if (d.size() != 3) {
    throw std::runtime_error("'x' must have 3 dimensions.");
  }
  
  if (d[0] != n || d[1] != n || d[2] != p)
    return 0;
  
  MatrixXd TrailC = MatrixXd::Map(Trail.begin(), nn, 3);
  
  MatrixXd NijT1 = MatrixXd::Zero(n, n);
  MatrixXd NijT2 = MatrixXd::Zero(n, n);
  
  MatrixXd V = MatrixXd::Zero(2*n - 1, 2*n - 1);
  MatrixXd W = MatrixXd::Zero(2*n - 1, 2*n - 1);
  MatrixXd S = MatrixXd::Zero(2*n - 1, 2*n - 1);
  MatrixXd O1 = MatrixXd::Zero(2*n - 1, 2*n - 1);
  MatrixXd O2 = MatrixXd::Zero(p, p);
  MatrixXd Eij = MatrixXd::Zero(n, n);
  MatrixXd Vn = MatrixXd::Zero(p, 2*n - 1);
  MatrixXd HQ = MatrixXd::Zero(p, p);
  MatrixXd HQi = MatrixXd::Zero(p, p);
  MatrixXd HQc = VectorXd::Zero(p);
  MatrixXd Sig = MatrixXd::Zero(p, p);
  MatrixXd Sigc = MatrixXd::Zero(p, 2*n - 1);
  
  VectorXd xk = VectorXd::Map(xe.begin(), 2*n + p - 1);
  Map<VectorXd> a(xk.data(), n);
  Map<VectorXd> b(xk.data() + n, n - 1);
  Map<VectorXd> c(xk.data() + 2*n - 1, p);
  
  VectorXd xkCI = VectorXd::Zero(2*n + p - 1);
  Map<VectorXd> abCI(xkCI.data(), 2*n - 1);
  Map<VectorXd> cCI(xkCI.data() + 2*n - 1, p);
  
  double diff(0.0), temp(0.0), v2nn(0.0);
  int i(0), j(0);
  
  static const double inv_sqrt_2pi = 0.3989422804014327;

  double bw1 = t - 10*h1 > 0 ? t - 10*h1 : 0;
  double bw2 = t + 10*h1 < 1 ? t + 10*h1 : 1;
  double bw3 = t - 10*h2 > 0 ? t - 10*h2 : 0;
  double bw4 = t + 10*h2 < 1 ? t + 10*h2 : 1;
  
  for (int z = 0; z < nn; z++) {
    if (TrailC(z, 2) > bw1 && TrailC(z, 2) < bw2) {
      i = TrailC(z, 0);
      j = TrailC(z, 1);
      diff = (TrailC(z, 2) - t) / h1;
      NijT1(i - 1, j - 1) += std::pow(inv_sqrt_2pi / h1 * std::exp(-0.5 * diff * diff), 2);
    }
    if (TrailC(z, 2) > bw3 && TrailC(z, 2) < bw4) {
      i = TrailC(z, 0);
      j = TrailC(z, 1);
      diff = (TrailC(z, 2) - t) / h2;
      NijT2(i - 1, j - 1) += std::pow(inv_sqrt_2pi / h2 * std::exp(-0.5 * diff * diff), 2);
    }
  }
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n - 1; j++) {
      if (i != j) {
        temp = a(i) + b(j);
        for (int k = 0; k < p; k++) {
          temp += c(k) * Zij[i + j*n + k*n*n];
        }
        Eij(i, j) = exp(temp);
      }
    }
    temp = a(i);
    for (int k = 0; k < p; k++) {
      temp += c(k) * Zij[i + (n-1)*n + k*n*n];
    }
    Eij(i, n - 1) = exp(temp);
  }
  
  
  for (int i = 0; i < n; i++) {
    temp = 0;
    for (int j = 0; j < n - 1; j++) {
      temp += Eij(i, j);
      V(i, n + j) = Eij(i, j) / (n - 1);
    }
    v2nn += Eij(i, n - 1) / (n - 1);
    V(i, i) = (temp + Eij(i, n - 1)) / (n - 1);
  }
  
  for (int j = 0; j < n - 1; j++) {
    temp = 0;
    for (int i = 0; i < n; i++) {
      temp += Eij(i, j);
      V(n + j, i) = Eij(i, j) / (n - 1);
    }
    V(n + j, n + j) = temp / (n - 1);
  }
  
  for (int i = 0; i < 2*n - 1; i++) {
    for (int j = 0; j < 2*n - 1; j++) {
      S(i, j) = -1 / v2nn;
    }
  }
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      S(i, j) = (i == j) / V(i, i) + 1/v2nn;
    }
  }
  
  for (int i = 0; i < n - 1; i++) {
    for (int j = 0; j < n - 1; j++) {
      S(i + n, j + n) = (i == j) / V(i + n, i + n) + 1/v2nn;
    }
  }
  
  for (int i = 0; i < n; i++) {
    temp = 0;
    for (int j = 0; j < n - 1; j++) {
      temp += NijT1(i, j);
      W(i, j + n) = NijT1(i, j);
    }
    temp += NijT1(i, n - 1);
    W(i, i) = temp;
  }
  
  for (int j = 0; j < n - 1; j++) {
    temp = 0;
    for (int i = 0; i < n; i++) {
      temp += NijT1(i, j);
      W(j + n, i) = NijT1(i, j);
    }
    W(j + n, j + n) = temp;
  }
  
  O1 = ((1000*S) * W * (1000*S)) / (n * n*1000000);
  
  for (int i = 0; i < 2*n - 1; i++) {
    abCI(i) = std::sqrt(O1(i, i));
  }
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (j != i) {
        for (int k = 0; k < p; k++)
          Vn(k, i) += Zij[i + j*n + k*n*n] * Eij(i, j);
      }
    }
  }
  
  for (int j = 0; j < n - 1; j++) {
    for (int i = 0; i < n; i++) {
      if (j != i) {
        for (int k = 0; k < p; k++)
          Vn(k, j + n) += Zij[i + j*n + k*n*n] * Eij(i, j);
      }
    }
  }
  
  Vn = Vn / n;
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i != j) {
        for (int k = 0; k < p; k++)
          HQc(k) = Zij[i + j*n + k*n*n];
        HQ += HQc * HQc.transpose() * Eij(i, j);
      }
    }
  }
  
  HQ = HQ / n + Vn * S * Vn.transpose();
  
  Sigc = Vn * S;
  
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < n; j++) {
      if (i != j) {
        for (int k = 0; k < p; k++)
          HQc(k) = Zij[i + j*n + k*n*n] - Sigc(k, i) - Sigc(k, j);
        Sig += HQc * HQc.transpose() * NijT2(i, j); 
      }
    }
  }

  Sig = Sig / n;
  
  HQi = HQ.partialPivLu().inverse();
  O2 = HQi * Sig * HQi;
  
  for (int i = 0; i < p; i++) {
    cCI(i) = std::sqrt(O2(i, i) / n);
  }
  
  return(Rcpp::List::create(Rcpp::_["O1"] = O1, Rcpp::_["S"] = S, Rcpp::_["W"] = W, Rcpp::_["V"] = V, Rcpp::_["vn"] = v2nn, Rcpp::_["xkCI"] = xkCI, Rcpp::_["Eij"] = Eij));
}
