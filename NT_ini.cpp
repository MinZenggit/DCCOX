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
NumericVector NewtonMC_ini(NumericMatrix Trail, NumericVector Zij, NumericVector xkIni, double t, double h1, double h2, int n, int nn, int p) {
    
    if (Rf_isNull(Zij.attr("dim"))) {
        throw std::runtime_error("'x' does not have 'dim' attibute.");
    }
    Rcpp::Dimension d = Zij.attr("dim");
    if (d.size() != 3) {
        throw std::runtime_error("'x' must have 3 dimensions.");
    }
    
    if (d[0] != n || d[1] != n || d[2] != p)
        return 0;
    
    MatrixXd NijT1 = MatrixXd::Zero(n, n);
    MatrixXd NijT2 = MatrixXd::Zero(n, n);
    VectorXd f1 = VectorXd::Zero(2*n + p - 1);
    VectorXd xk = VectorXd::Map(xkIni.begin(), 2*n + p - 1);
    Map<VectorXd> a(xk.data(), n);
    Map<VectorXd> b(xk.data() + n, n - 1);
    Map<VectorXd> c(xk.data() + 2*n - 1, p);
    
    VectorXd dd = VectorXd::Zero(2*n + p - 1);
    VectorXd f = VectorXd::Zero(2*n + p - 1);
    MatrixXd D = MatrixXd::Zero(2*n + p - 1, 2*n + p - 1);
    MatrixXd Eij = MatrixXd::Zero(n, n);
    
    MatrixXd TrailC = MatrixXd::Map(Trail.begin(), nn, 3);
    
    int i(0), j(0);
    double conv(1.0), mf1(10.0), temp(0.0);
    double bw1 = t - 10*h1 > 0 ? t - 10*h1 : 0;
    double bw2 = t + 10*h1 < 1 ? t + 10*h1 : 1;
    double bw3 = t - 10*h2 > 0 ? t - 10*h2 : 0;
    double bw4 = t + 10*h2 < 1 ? t + 10*h2 : 1;
    
    static const double inv_sqrt_2pi = 0.3989422804014327;
    double diff(0.0);
    
    for (int z = 0; z < nn; z++) {
        i = TrailC(z, 0);
        j = TrailC(z, 1);
        NijT1(i - 1, j - 1) += 1;
        NijT2(i - 1, j - 1) += 1;
    }
    
    
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++){
            if (j != i) {
                f1(i) += NijT1(i, j);
            }
        }
    }
    
    for (int j = 0; j < n - 1; j++) {
        for (int i = 0; i < n; i++) {
            if (i != j) {
                f1(n + j) += NijT1(i, j);
            }
        }
    }
    
    for (int k = 0; k < p; k++) {
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                if (i != j) {
                    f1(2*n - 1 + k) += Zij[i + j*n + k*n*n] * NijT2(i, j);
                }
            }
        }
    }

    
    for (int i = 0; i < 2*n - 1; i++) {
        if (f1(i) != 0 && mf1 > f1(i))
            mf1 = f1(i);
    }
    
    for (int i = 0; i < 2*n - 1; i++) {
        if (f1(i) == 0)
            f1(i) = mf1;
    }
    
    while (conv > 0.0001) {
        
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
            f(i) = 0;
            for (int j = 0; j < n; j++) {
                if (j != i) {
                    f(i) -= Eij(i, j);
                }
            }
            f(i) += f1(i);
        }
        
        for (int j = 0; j < n - 1; j++) {
            f(n + j) = 0;
            for (int i = 0; i < n; i++) {
                if (i != j) {
                    f(n + j) -= Eij(i, j);
                }
            }
            f(n + j) += f1(n + j);
        }
        
        for (int k = 0; k < p; k++) {
            f(2*n - 1 + k) = 0;
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        f(2*n - 1 + k) -= Zij[i + j*n + k*n*n] * Eij(i, j);
                    }
                }
            }
            f(2*n - 1 + k) += f1(2*n - 1 + k);
        }

        
        for (int i = 0; i < n; i++) {
            D(i, i) = 0;
            for (int k = 0; k < p; k++)
                D(i, 2*n - 1 + k) = 0;
            for (int j = 0; j < n; j++) {
                if (i != j) {
                    D(i, i) += Eij(i, j);
                    for (int k = 0; k < p; k++)
                        D(i, 2*n - 1 + k) += Zij[i + j*n + k*n*n] * Eij(i, j);
                }
            }
            for (int k = 0; k < p; k++)
                D(2*n - 1 + k, i) = D(i, 2*n - 1 + k);
        }
        
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n - 1; j++) {
                if (i != j) {
                    D(i, j + n) = Eij(i, j);
                }
            }
        }
        
        for (int j = 0; j < n - 1; j++) {
            for (int i = 0; i < n; i++) {
                if (i != j) {
                    D(j + n, i) = Eij(i, j);
                }
            }
        }
        
        for (int j = 0; j < n - 1; j++) {
            D(j + n, j + n) = 0;
            for (int k = 0; k < p; k++)
                D(j + n, 2*n - 1 + k) = 0;
            for (int i = 0; i < n; i++) {
                if (i != j) {
                    D(j + n, j + n) += Eij(i, j);
                    for (int k = 0; k < p; k++)
                        D(j + n, 2*n - 1 + k) += Zij[i + j*n + k*n*n] * Eij(i, j);
                }
            }
            for (int k = 0; k < p; k++)
                D(2*n - 1 + k, j + n) = D(j + n, 2*n - 1 + k);
        }
        
        for (int k1 = 0; k1 < p; k1++) {
            for (int k2 = 0; k2 < p; k2++)
                D(2*n - 1 + k1, 2*n - 1 + k2) = 0;
        }
        
        for (int k1 = 0; k1 < p; k1++) {
            for (int k2 = 0; k2 < p; k2++) {
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        if (i != j) {
                            D(2*n - 1 + k1, 2*n - 1 + k2) += Zij[i + j*n + k1*n*n] * Zij[i + j*n + k2*n*n] * Eij(i, j);
                        }
                    }
                }
            }
        }
        
        dd = D.partialPivLu().solve(f);
        xk += 0.1*dd; // NEED JUDGE
        conv = dd.norm();
        Rcpp::Rcout << conv << "\n";
    }
     
    
    return Rcpp::wrap(xk);
    // return Rcpp::List::create(Rcpp::Named("xk") = xk , Rcpp::Named("D") = D, Rcpp::Named("f") = f, Rcpp::Named("Eij") = Eij);
}