// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
#include <cmath>
using namespace Rcpp;
#define EIGEN_NO_DEBUG
#define EIGEN_DONT_PARALLELIZE

using Rcpp::NumericMatrix;
using Rcpp::NumericVector;
using Eigen::Map;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// [[Rcpp::export]]
NumericVector NewtonMCHomo(NumericMatrix Trail, NumericVector Zij, double t, double h2, int n, int nn, int p) {
    
    if (Rf_isNull(Zij.attr("dim"))) {
        throw std::runtime_error("'x' does not have 'dim' attibute.");
    }
    Rcpp::Dimension d = Zij.attr("dim");
    if (d.size() != 3) {
        throw std::runtime_error("'x' must have 3 dimensions.");
    }
    
    if (d[0] != n || d[1] != n || d[2] != p)
        return 0;
    
    MatrixXd NijT2 = MatrixXd::Zero(n, n);
    VectorXd f1 = VectorXd::Zero(p);
    VectorXd xk = VectorXd::Zero(p);
    Map<VectorXd> c(xk.data(), p);
    
    VectorXd dd = VectorXd::Zero(p);
    VectorXd f = VectorXd::Zero(p);
    MatrixXd D = MatrixXd::Zero(p, p);
    MatrixXd Eij = MatrixXd::Zero(n, n);
    
    MatrixXd TrailC = MatrixXd::Map(Trail.begin(), nn, 3);
    
    int i(0), j(0);
    double conv(1.0), mf1(10.0), temp(0.0);
    double bw3 = t - 10*h2 > 0 ? t - 10*h2 : 0;
    double bw4 = t + 10*h2 < 1 ? t + 10*h2 : 1;
    
    static const double inv_sqrt_2pi = 0.3989422804014327;
    double diff(0.0);
    // 计算 \int_{0}^1 \mathcal K_{h1}(s-t) dN_{ij}(s)
    for (int z = 0; z < nn; z++) {
        if (TrailC(z, 2) > bw3 && TrailC(z, 2) < bw4) {
            i = TrailC(z, 0);
            j = TrailC(z, 1);
            diff = (TrailC(z, 2) - t) / h2;
            NijT2(i - 1, j - 1) += inv_sqrt_2pi / h2 * std::exp(-0.5 * diff * diff);
        }
    }
    // 对于z的每一维度 k，计算 \sum_{i=1}^n \sum_{j \neq i} \int_{0}^1 Z_{ij}(s) \mathcal K_{h_2}(s-t) d N_{ij}(s)
    for (int k = 0; k < p; k++) {
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                if (i != j) {
                    f1(k) += Zij[i + j*n + k*n*n] * NijT2(i, j);
                }
            }
        }
    }

    while (conv > 0.0001) {
        
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    temp = 0;
                    for (int k = 0; k < p; k++) {
                        temp += c(k) * Zij[i + j*n + k*n*n];
                    }
                    Eij(i, j) = exp(temp);
                    // exp(Z_ij*c)
                }
            }
        }
        
        for (int k = 0; k < p; k++) {
            f(k) = 0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        f(k) -= Zij[i + j*n + k*n*n] * Eij(i, j);
                    }
                }
            }
            f(k) += f1(k);
        }
        for (int k1 = 0; k1 < p; k1++) {
            for (int k2 = 0; k2 < p; k2++) {
                D(k1, k2) = 0;
              //重要的；
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        if (i != j) {
                            D(k1, k2) += Zij[i + j*n + k1*n*n] * Zij[i + j*n + k2*n*n] * Eij(i, j);
                        }
                    }
                }
            }
        }
        
        dd = D.partialPivLu().solve(f);
        // Rcpp::Rcout << dd.norm() << std::endl;
        xk += 0.05*dd;
        conv = dd.norm();
    }
     
    
    return Rcpp::wrap(xk);
}