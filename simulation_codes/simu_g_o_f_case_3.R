rm(list = ls())
#  simulation for goodness of fit, case 1: alpha_i = beta_i = f.
library(matrixcalc)
library(dplyr)
library(Rcpp)
library(RcppEigen)
library(ggplot2)
# Sys.setenv("PKG_CPPFLAGS" = "-march=native")
sourceCpp(file = "cpp/NT.cpp")
sourceCpp(file = "cpp/NT_homo.cpp")
sourceCpp(file = "cpp/SimSet_c.cpp", rebuild = T)
sourceCpp(file = "cpp/CI.cpp")
# source(file = "NewtonMCHomo_R.R")
set.seed(100)
n = 100
p_true = 3 # covariates dimention
# 4 groups
t = 0
b = 1

fg <- function(t,i){
  sin(2*pi*t)/3
}


xkkk = xkkk2 = list()
cont = 1
nnn = 0


# zij_true <- array(0, c(n, n, p_true))
# random_data <- rbinom(n * n, 1, 0.5)  # 产生一些小的噪声
# zij_true[,,1] <- matrix(random_data, n, n)
# correlation_factor <- rnorm(n * n, mean = 0, sd = 0.1)
# zij_true[,,2] <- zij_true[,,1] + matrix(correlation_factor, n, n)
# zij_true[,,2] <- ifelse(zij_true[,,2] > 1, 1, ifelse(zij_true[,,2] < 0, 0, zij_true[,,2]))
# correlation_factor <- rnorm(n * n, mean = 0, sd = 0.1)  # 产生一些小的噪声
# zij_true[,,3] <- zij_true[,,2] + matrix(correlation_factor, n, n)
# zij_true[,,3] <- ifelse(zij_true[,,3] > 1, 1, ifelse(zij_true[,,3] < 0, 0, zij_true[,,2]))
zij_true = array(rbinom(n*n*p_true, 1, 0.5), c(n,n,p_true))

trail_sim = SimSetC(n, b, zij_true)
trail_sim = as.data.frame(trail_sim)
rv = order(trail_sim[,3])
trail_sim = trail_sim[rv,]
nn = nrow(trail_sim)
colnames(trail_sim) <- c("s", "r", "t")

# p=1
# zij = array(0, c(n,n,p))
# zij[,,1] <- zij_true[,,1]
p=3
zij <- zij_true
# zij = array(0, c(n, n, p))
# zij[1:4, 1:(n/3), ] = 1
ones = array(1, c(n, n, 1))
zij_new = abind::abind(zij, ones, along = 3)
#simulate data with reject and accept method; though constant z_ij


# count the total events
Nij = matrix(0, nrow = n, ncol = n)
for (z in 1:nn) {
  i = trail_sim[z, 1]
  j = trail_sim[z, 2]
  Nij[i, j] = Nij[i, j] + 1
}

colMeans(Nij[1:50, ])
colMeans(Nij[51:100, ])
# Newton Method -----------------------------------------------------------


xkk = matrix(rep(0, 2*n + p - 1)) # for all parameters
xkk2 = xkk3 = matrix(rep(0, p+1))
for (t in seq(0.1,0.9,0.05)) {
  
  print(t)
  h1 = 0.3*n^(-0.1)
  h2 = 0.015*n^(-0.2)
  
  xk1 = NewtonMC(as.matrix(trail_sim), zij, t, h1, h2, n, nn, p)
  xk2 = NewtonMCHomo(as.matrix(trail_sim), zij_new, t, h2, n, nn, p+1)
  # xk3 = NewtonMCHomo_R(as.matrix(trail_sim), zij_new, t, h2, n, nn, p+1)
  
  xkk = cbind(xkk, matrix(xk1))
  xkk2 = cbind(xkk2, matrix(xk2))
  # xkk3 = cbind(xkk3, matrix(xk3))
}

xkk = xkk[, -1]
xkkHomo = xkk2[, -1]
# xkkHomo2 = xkk3[, -1]


pk=1
p1 = data.frame(t = rep(seq(0.1,0.9,0.05),2),
                y = c(xkk[2*n-1+pk,], xkkHomo[pk,]),
                het = c(rep("y", 17), rep("n", 17)))

ggplot(p1, aes(x = t, y = y, group = het, color = het, fill = het)) +
  geom_line(size = 0.7) +
  stat_function(fun = fg, args = list(i = pk), color = "black", size = 0.7) +
  scale_color_manual(values=c("#619CFF", "red")) +
  xlab(expression(italic("t"))) +
  ylab(expression(hat(italic(gamma)))) +
  theme_bw() +
  theme(legend.position = "none") -> p1
plot(p1)



########
# Ours
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
int_exp_ij <- exp_ij <- zij_list <- gs_list <- list()
for (i in 1:n) {
  for (j in 1:n) {
    exp_ij[[paste0(i, "->", j)]] <- local({
      current_i <- i
      current_j <- j
      current_n <- n
      as_f <- a_b_g_hat[[current_i]]
      bs_f <- a_b_g_hat[[current_n + current_j]]
      for(pp in 1:p){
        zij_list[[pp]] <- zij[current_i, current_j, pp]
        gs_list[[pp]] <- a_b_g_hat[[2 * current_n - 1 +pp]]
      }
      
      function(s) {
        if (current_j != current_n) {
          temp = as_f(s) + bs_f(s)
          for(pp in 1:p){
            temp = temp + zij_list[[pp]] * gs_list[[pp]](s)
          }
          return(exp(temp))
        }
        if (current_j == current_n) {
          temp = as_f(s) + 0
          for(pp in 1:p){
            temp = temp + zij_list[[pp]] * gs_list[[pp]](s)
          }
          return(exp(temp))
        }
      }
    })
    int_exp_ij[[paste0(i, "->",j)]] <- local({
      current_i = i
      current_j = j
      exp_ij_f <- exp_ij[[paste0(current_i, "->",current_j)]]
      function(t){
        integrate(exp_ij_f, lower = 0.1, upper = t)$value
      }
    })
  }
}



# iseq <- c(1, 25, 51, 75)
iseq <- c(1:n)
es_Ni_list <- list()
for( i in iseq){
  print(i)
  trail_sim %>% filter(s == i, t>=0.1, t<0.9) -> dd
  dd[order(dd$t), ] -> dd
  dd$N = 1:nrow(dd)
  if(nrow(dd) > 100){
    ind = sample(1:nrow(dd), 100, replace = F) %>% sort()
  }else{
    ind = 1:nrow(dd)
  }
  es_Ni <- c()
  for(t in dd$t[ind]){
    int_exp_i = 0
    for(j in 1:n){
      if(j != i)
        int_exp_i = int_exp_i + int_exp_ij[[paste0(i, "->", j)]](t)
    }
    es_Ni <- c(es_Ni, int_exp_i)
  }
  es_Ni_list[[i]] = rbind(N = dd$N[ind], es_Ni = es_Ni)
}
save(es_Ni_list, file = "simu_results/goodness_of_fit/es_Ni_case_3.rdata")
re_i <- c()
for(i in iseq){
  re <- cbind(t(es_Ni_list[[i]]), i)
  colnames(re) <- c("Number of events of senders", "Esimtated cumulative intensity", "Sender")
  re_i <- rbind(re_i, re)
}
re_i %>% as.data.frame() -> re_i
ggplot(re_i, aes(x = `Number of events of senders`, y = `Esimtated cumulative intensity`, color = `Sender`)) +
  geom_point(size = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  labs(title = "Our's methods: sender",
       x = "Number of events of senders",
       y = "Estimated cumulative intensity") +
  theme_minimal() -> p1
plot(p1)
#
#
#
# plot(seq(0.1, 0.9, 0.05), sapply(seq(0.1,0.9,0.05), fs, i=1))
# points(seq(0.1, 0.9, 0.05), xkk[2*n, ], col = "blue")
# points(seq(0.1, 0.9, 0.05), xkk[2*n+2, ], col = "green")

#
es_Nj_list <- list()
for( j in iseq){
  cat(j, "\n")
  trail_sim %>% filter(r == j, t>=0.1, t<0.9) -> dd
  dd[order(dd$t), ] -> dd
  dd$N = 1:nrow(dd)
  if(nrow(dd) > 100){
    ind = sample(1:nrow(dd), 100, replace = F) %>% sort()
  }else{
    ind = 1:nrow(dd)
  }
  es_Nj <- c()
  for(t in dd$t[ind]){
    int_exp_j = 0
    for(i in 1:n){
      if(i != j )
        int_exp_j = int_exp_j + int_exp_ij[[paste0(i, "->", j)]](t)
    }
    es_Nj <- c(es_Nj, int_exp_j)
  }
  es_Nj_list[[j]] = rbind(N = dd$N[ind], es_Nj = es_Nj)
}
save(es_Nj_list, file = "simu_results/goodness_of_fit/es_Nj_case_3.rdata")

##################
# Homo
homo_hat <- list()
for (pk in 1:(p + 1)) {
  homo_hat[[pk]] <- local({
    current_pk <- pk  # 捕获当前循环的值
    function(s) {
      fit <- smooth.spline(seq(0.1, 0.9, 0.05), xkkHomo[current_pk, ])
      return(predict(fit, s)$y)
    }
  })
}
# homo_hat[[1]] <- homo_hat[[2]]  <-  function(s){return(as.numeric(s-0.5))}

int_exp_ij_homo <- exp_ij_homo <- zij_list <- gs_list <- list()
for (i in 1:n) {
  for (j in 1:n) {
    exp_ij_homo[[paste0(i, "->", j)]] <- local({
      current_i <- i
      current_j <- j
      current_n <- n
      for(pp in 1:(p+1)){
        zij_list[[pp]] <- zij_new[current_i, current_j, pp]
        gs_list[[pp]] <- homo_hat[[pp]]
      }
      
      function(s) {
        if (current_j != current_n) {
          temp = 0
          for(pp in 1:(p+1)){
            temp = temp + zij_list[[pp]] * gs_list[[pp]](s)
          }
          return(exp(temp))
        }
        if (current_j == current_n) {
          temp = 0
          for(pp in 1:(p+1)){
            temp = temp + zij_list[[pp]] * gs_list[[pp]](s)
          }
          return(exp(temp))
        }
      }
    })
    int_exp_ij_homo[[paste0(i, "->",j)]] <- local({
      current_i = i
      current_j = j
      exp_ij_f <- exp_ij_homo[[paste0(current_i, "->",current_j)]]
      function(t){
        integrate(exp_ij_f, lower = 0.1, upper = t)$value
      }
    })
  }
}

es_Ni_homo_list <- list()

for( i in iseq){
  print(i)
  trail_sim %>% filter(s == i, t>=0.1, t<0.9) -> dd
  dd[order(dd$t), ] -> dd
  dd$N = 1:nrow(dd)
  if(nrow(dd) > 100){
    ind = sample(1:nrow(dd), 100, replace = F) %>% sort()
  }else{
    ind = 1:nrow(dd)
  }
  es_Ni_homo <- c()
  for(t in dd$t[ind]){
    int_exp_i_homo = 0
    for(j in 1:n){
      if(j !=i )
        int_exp_i_homo = int_exp_i_homo + int_exp_ij_homo[[paste0(i, "->", j)]](t)
    }
    es_Ni_homo <- c(es_Ni_homo, int_exp_i_homo)
  }
  es_Ni_homo_list[[i]] = rbind(N = dd$N[ind], es_Ni_homo = es_Ni_homo)
}
save(es_Ni_homo_list, file = "simu_results/goodness_of_fit/homo_es_Ni_case_3.rdata")
# 
re_i <- c()
for(i in iseq){
  re <- cbind(t(es_Ni_homo_list[[i]]), i)
  colnames(re) <- c("Number of events of senders", "Esimtated cumulative intensity", "Sender")
  re_i <- rbind(re_i, re)
}
re_i %>% as.data.frame() -> re_i
ggplot(re_i, aes(x = `Number of events of senders`, y = `Esimtated cumulative intensity`, color = `Sender`)) +
  geom_point(size = 0.5) +
  geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
  labs(title = "K's methods: sender",
       x = "Number of events of senders",
       y = "Estimated cumulative intensity") +
  theme_minimal() -> p1
plot(p1)
# plot(seq(0.1, 0.9, 0.05), sapply(seq(0.1,0.9,0.05), fg, i=1))
# points(seq(0.1, 0.9, 0.05), xkkHomo[1, ], col = "blue")
# points(seq(0.1, 0.9, 0.05), xkkHomo[2, ], col = "red")



es_Nj_homo_list <- list()
for( j in iseq){
  cat(j, "\n")
  trail_sim %>% filter(r == j, t>=0.1, t<0.9) -> dd
  dd[order(dd$t), ] -> dd
  dd$N = 1:nrow(dd)
  if(nrow(dd) > 100){
    ind = sample(1:nrow(dd), 100, replace = F) %>% sort()
  }else{
    ind = 1:nrow(dd)
  }
  es_Nj_homo <- c()
  for(t in dd$t[ind]){
    int_exp_j_homo = 0
    for(i in 1:n){
      if(i != j )
        int_exp_j_homo = int_exp_j_homo + int_exp_ij_homo[[paste0(i, "->", j)]](t)
    }
    es_Nj_homo <- c(es_Nj_homo, int_exp_j_homo)
  }
  es_Nj_homo_list[[j]] = rbind(N = dd$N[ind], es_Nj_homo = es_Nj_homo)
}
save(es_Nj_homo_list, file = "simu_results/goodness_of_fit/homo_es_Nj_case_3.rdata")
