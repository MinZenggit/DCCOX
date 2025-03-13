
# rm(list = ls())
set.seed(100)
# n=100 seed 100
p = 2 # covariates dimention
n = 100
K = 5
cont = 1

Total_cont <- 1

seq(0.05, 0.5, by = 0.05) -> h1_set_c
c(seq(0.002, 0.03, by = 0.004), 0.04, 0.12) -> h2_set_c

h1_set <- h1_set_c*n^(-0.1)
h2_set <- h2_set_c*n^(-0.2)


library(matrixcalc)
library(Rcpp)
library(RcppEigen)
sourceCpp(file = "cpp/NT_cv.cpp", verbose = TRUE, rebuild = TRUE)
sourceCpp(file = "cpp/SimSet.cpp", verbose = TRUE, rebuild = TRUE)

# 4 groups
t = 0
cs = 0.5
shift = cs*log(n)

# s for sender; r for receiver
fs <- function(t, i) {
  if (i < n/2) return(2.5-shift+sin(2*pi*t))
  return(1.5+t/2-shift)
} 

fr <- function(t, i) {
  if (i < n/2) return(2.5-shift+cos(2*pi*t))
  if (i < n) return(1.5+t/2-shift)
  return(0)
} 

fg <- function(t, i) {
  return(sin(2*pi*t)/3)
} 

X_true1 = matrix(sapply(seq(0.1,0.9,0.05), fs, i = 1), nrow = (n/2-1), ncol = 17, byrow = TRUE)
X_true2 = matrix(sapply(seq(0.1,0.9,0.05), fs, i = n/2), nrow = (n/2+1), ncol = 17, byrow = TRUE)
X_true3 = matrix(sapply(seq(0.1,0.9,0.05), fr, i = 1), nrow = (n/2-1), ncol = 17, byrow = TRUE)
X_true4 = matrix(sapply(seq(0.1,0.9,0.05), fs, i = n/2), nrow = (n/2), ncol = 17, byrow = TRUE)
X_true5 = matrix(sapply(seq(0.1,0.9,0.05), fg, i = 1), nrow = p, ncol = 17, byrow = TRUE)
X_true = rbind(X_true1, X_true2, X_true3, X_true4, X_true5)

PE_m <- matrix(0, nrow = length(h1_set), ncol = length(h2_set))
MISEs <- array(0, dim = c(length(h1_set), length(h2_set), 6))
while(cont <= Total_cont) {
  cat("\n cont =", cont, "\ ")
  # Generate data
  zij = array(rnorm(n*n*p), c(n,n,p))
  
  #simulate data with reject and accept method; though constant z_ij
  trail_sim = SimSetC(n, 0.5, zij)
  trail_sim = as.data.frame(trail_sim)
  
  rv = order(trail_sim[,3])
  trail_sim = trail_sim[rv,]
  nn = nrow(trail_sim)
  colnames(trail_sim) <- c("s", "r", "t")
  
  # count the total events
  Nij = matrix(0, nrow = n, ncol = n)
  for (z in 1:nn) {
    i = trail_sim[z, 1]
    j = trail_sim[z, 2]
    Nij[i, j] = Nij[i, j] + 1
  }
  
  # 生成CV矩阵
  CVij <- array(1, dim = c(n, n, K))
  KK <- c()
  for(ip in 1:(n/K)){
    tt <- sample(1:n, n, replace = F)
    Subgroup_K <- array(1, dim = c(K, n))
    for(ik in 1:K){
      for(il in 1:K){
        b=((ik-1 + il-1)*(n/K)+1)
        e=((ik+il-1)*(n/K))
        if(b > n){
          b = b-n
        }
        if(e > n){
          e = e-n
        }
        Subgroup_K[ik, tt[b:e]] = il
      }
    }
    KK <- rbind(KK, Subgroup_K)
  }
  
  for(i in 1: n){
    for(k in 1: K){
      CVij[i, which(KK[i, ]==k), k] = 0
    }
  }
  
  
  # For different h
  for(ii in 1: length(h1_set)){
   for(jj in 1: length(h2_set)){
      h1 <- h1_set[ii]
      h2 <- h2_set[jj]
      cat("\n h1 =", h1, ", h2 =", h2, ".\n k = ")
      
      
      # For fixed h, do CV
      # splines for the PE calculation
      PE <- numeric(K)
      for(k in 1: K){
        # for each k，eatimete a, b, g in seq points
        cat(k, " ,")
        xkk = matrix(rep(0, 2*n + p - 1)) # for all parameters
        for (t in seq(0.1,0.9,0.05)) {
          xk1 = NewtonMC_CV(as.matrix(trail_sim), zij, CVij[,,k],t, h1, h2, n, nn, p)
          xkk = cbind(xkk, matrix(xk1))
        }
        xkk = xkk[, -1]
        # a, b, g: spline
        a_b_g_hat <- list()
        for (pk in 1:(2 * n + p - 1)) {
          a_b_g_hat[[pk]] <- local({
            current_pk <- pk  # 捕获当前循环的值
            function(s) {
              fit <- smooth.spline(seq(0.1, 0.9, 0.05), xkk[current_pk, ])
              return(predict(fit, s)$y)
            }
          })
        }
        # 计算积分，对于每一条边，注意beta_n始终为0
        int_exp_ij <- exp_ij <- list()
        for (i in 1:n) {
          for (j in 1:n) {
            if(CVij[i, j, k] == 0){
              exp_ij[[paste0(i, "->", j)]] <- local({
                current_i <- i
                current_j <- j
                current_n <- n
                zij_1 <- zij[current_i, current_j, 1]
                zij_2 <- zij[current_i, current_j, 2]
                as_f <- a_b_g_hat[[current_i]]
                bs_f <- a_b_g_hat[[n + current_j]]
                gs_1f <- a_b_g_hat[[2 * current_n]]
                gs_2f <- a_b_g_hat[[2 * current_n + 1]]
                
                function(s) {
                  if (current_j != current_n) {
                    return(exp(as_f(s) + bs_f(s) +
                                 zij_1 * gs_1f(s) + 
                                 zij_2 * gs_2f(s)))
                  }
                  if (current_j == current_n) {
                    return(exp(as_f(s) + 0 +
                                 zij_1 * gs_1f(s) + 
                                 zij_2 * gs_2f(s)))
                  }
                }
              })
              int_exp_ij[[paste0(i, "->",j)]] <- local({
                current_i = i
                current_j = j
                exp_ij_f <- exp_ij[[paste0(current_i, "->",current_j)]]
                function(t){
                  integrate(exp_ij_f, lower = 0, upper = t)$value
                }
              })
            }
          }
        }
        for(e in 1:nrow(trail_sim)){
          s <- trail_sim$s[e]
          r <- trail_sim$r[e]
          t <- trail_sim$t[e]
          if(CVij[s, r, k] == 0){
            PE[k] <- PE[k] + (1 - int_exp_ij[[paste0(s, "->",r)]](t))^2
          }
        }
      }
      PE_m[ii, jj] = PE_m[ii, jj] + sum(PE)/Total_cont
      
      #MISEs of a_1
      MISEs[ii, jj, 1] <- MISEs[ii, jj, 1] + integrate(function(s){
        (a_b_g_hat[[1]](s) - fs(s, 1))^2}, 0, 1)$value/Total_cont
      #MISEs of a_n/2+1
      MISEs[ii, jj, 2] <- MISEs[ii, jj, 2] + integrate(function(s){
        (a_b_g_hat[[n/2+1]](s) - fs(s, n/2+1))^2}, 0, 1)$value/Total_cont
      
      #MISEs of b_1
      MISEs[ii, jj, 3] <- MISEs[ii, jj, 3] + integrate(function(s){
        (a_b_g_hat[[n/2+1]](s) - fr(s, 1))^2}, 0, 1)$value/Total_cont
      #MISEs of b_n/2+1
      MISEs[ii, jj, 4] <- MISEs[ii, jj, 4] + integrate(function(s){
        (a_b_g_hat[[n/2+1]](s) - fr(s, n/2+1))^2}, 0, 1)$value/Total_cont
      
      #MISEs of g_1
      MISEs[ii, jj, 5] <- MISEs[ii, jj, 5] + integrate(function(s){
        (a_b_g_hat[[2*n]](s) - fg(s))^2}, 0, 1)$value/Total_cont
      #MISEs of g_2
      MISEs[ii, jj, 6] <- MISEs[ii, jj, 6] + integrate(function(s){
        (a_b_g_hat[[2*n+1]](s) - fg(s))^2}, 0, 1)$value/Total_cont
      
     }
  }
  cont = cont + 1
}
rownames(PE_m) =paste0("h1=", round(h1_set_c, 3)) 
colnames(PE_m) =paste0("h2=", round(h2_set_c, 33)) 


my_colors <- colorRampPalette(c("blue", "red"))(256)
heatmap(PE_m, main="PE with different h1 and h2", col=my_colors, scale="column", Rowv=NA, Colv=NA)
# print(MISEs)
write.csv(PE_m, file = paste0("simu_results/simu_CV/n=", n, ".csv"))


