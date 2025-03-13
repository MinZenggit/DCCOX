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
NumericVector NewtonMC(NumericMatrix Trail, NumericVector Zij, double t, double h1, double h2, int n, int nn, int p) {
    
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
    VectorXd xk = VectorXd::Zero(2*n + p - 1);
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
        if (TrailC(z, 2) > bw1 && TrailC(z, 2) < bw2) {
            i = TrailC(z, 0);
            j = TrailC(z, 1);
            diff = (TrailC(z, 2) - t) / h1;
            //这个是每个边 event 关于时间t 做kernel加权的，用的是高斯kernel，
            // 计算 \int_{0}^1 \mathcal K_{h1}(s-t) dN_{ij}(s)
            NijT1(i - 1, j - 1) += inv_sqrt_2pi / h1 * std::exp(-0.5 * diff * diff);
        }
        // 除了h1，h2的不同，别的都一样
        if (TrailC(z, 2) > bw3 && TrailC(z, 2) < bw4) {
            i = TrailC(z, 0);
            j = TrailC(z, 1);
            diff = (TrailC(z, 2) - t) / h2;
            NijT2(i - 1, j - 1) += inv_sqrt_2pi / h2 * std::exp(-0.5 * diff * diff);
        }
    }
    
    // 对于每一个i，计算 \sum_{j \neq i} \int_{0}^1 \mathcal K_{h1}(s-t) dN_{ij}(s).
    // 对应的是算法中a_i(t)的分子部分
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++){
            if (j != i) {
                f1(i) += NijT1(i, j); // kernel加权平均结果，对于每一个边
            }
        }
    }
    // 对于每一个j，计算 \sum_{i \neq j} \int_{0}^1 \mathcal K_{h2}(s-t) dN_{ij}(s).
    // 对应的是算法中b_i(t)的分子部分
    for (int j = 0; j < n - 1; j++) {
        for (int i = 0; i < n; i++) {
            if (i != j) {
                f1(n + j) += NijT1(i, j);
            }
        }
    }
    // 对于z的每一维度 k，计算 \sum_{i=1}^n \sum_{j \neq i} \int_{0}^1 Z_{ij}(s) \mathcal K_{h_2}(s-t) d N_{ij}(s)
    // 对应算法中gamma 估计的前半部分。
    for (int k = 0; k < p; k++) {
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                if (i != j) {
                    f1(2*n - 1 + k) += Zij[i + j*n + k*n*n] * NijT2(i, j);
                }
            }
        }
    }
    
    //这一步是将==0的用最小正数替代
    for (int i = 0; i < 2*n - 1; i++) {
        if (f1(i) != 0 && mf1 > f1(i))
            mf1 = f1(i);
    }
    for (int i = 0; i < 2*n - 1; i++) {
        if (f1(i) == 0)
            f1(i) = mf1;
    }
    
  //上面这些计算是不需要迭代的。
  
  //下面的计算是需要每次迭代中进行更新的。  
    while (conv > 0.0001) {
      // 计算了 F(a_k, b_k, c_k)； 记做这里的f
      // 利用牛顿迭代法 F(a_k+1, b_k+1, c_k+1) = F(a_k, b_k, c_k) + dd
      // dd = - Jacobi(a_k, b_k, c_k)^{-1} * F(a_k, b_k, c_k)
      // 这里的Jacobi就是下面的D，然后dd是用求解线性方程的形式得到的。
        
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n - 1; j++) {
                if (i != j) {
                    temp = a(i) + b(j);
                    for (int k = 0; k < p; k++) {
                        temp += c(k) * Zij[i + j*n + k*n*n];
                    }
                    Eij(i, j) = exp(temp);
                    // exp(a_i+b_j+Z_ij*c)
                }
            }
            temp = a(i);
            for (int k = 0; k < p; k++) {
                temp += c(k) * Zij[i + (n-1)*n + k*n*n];
            }
            Eij(i, n - 1) = exp(temp);
            // exp(a_i + z_ij * c) // for identifiablity, set b_n = 0
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
        xk += dd;
        conv = dd.norm();
    }
     
    return Rcpp::wrap(xk);
}