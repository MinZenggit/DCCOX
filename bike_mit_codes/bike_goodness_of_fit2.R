library(RcppEigen)
library(Rcpp)
library(matrixcalc)
library(dplyr)
rm(list = ls())
setwd("/home/szj/DCCOX")
# setwd("C:/Users/micha/Downloads/dynamic_network")
sourceCpp(file = "cpp/NT_test.cpp")
sourceCpp(file = "cpp/NT_ini.cpp")
sourceCpp(file = "cpp/CI.cpp")
sourceCpp(file = "cpp/CI_alt.cpp")
sourceCpp(file = "cpp/NT_homo.cpp")

load(file = "bikedata/bike_data3.rdata")

# zij <- bike_data$zij[,,c(3,4)]
# trail <- bike_data$trail

nodes <- 1:100
zij <- bike_data$zij[nodes,nodes,c(1:4)]
zij[is.na(zij)] = 0
trail <- bike_data$trail %>% filter(x %in% nodes, y %in% nodes)
n = dim(zij)[1]
nn = nrow(trail)
p = dim(zij)[3]

h1 = 0.15*n^(-0.1)#0.04375991
h2 = 0.4*n^(-0.2)#0.04787323
tseq = seq(0.1,0.9,0.01)
# tseq = seq(0.05,0.95,0.01)[c(13, 17)]

xkk = matrix(rep(0, 2*n+p-1))

cat(h1, h2, "\n")
# xk1 = NewtonMC_ini(as.matrix(trail), zij, rep(0, 2*n + p - 1), 0.1, h1, h2, n, nn, p)

for (t in tseq) {
  xk2 = NewtonMC_test(as.matrix(trail), zij, rep(0.1, 2*n + p - 1), t, h1, h2, n, nn, p)
  xkk = cbind(xkk, matrix(xk2))
  cat(t, "\n")
}


xkk = xkk[, -1]
# xkkHomo = matrix(rep(0, p))
xkkHomo = matrix(rep(0, p+1))
ones = array(1, c(n, n, 1))
zij_ones = abind::abind(zij, ones, along = 3)


for (t in tseq) {
  # xk2 = NewtonMC_test(as.matrix(trail), zij, xk1, t, h1, h2, n, nn, 16)
  # xk3 = NewtonMCHomo(as.matrix(trail), zij_ones[,,5:6], t, h2, n, nn, p+1)
  # xk3 = NewtonMCHomo(as.matrix(trail), zij_ones[,,c(1,6)], t, h2, n, nn, 5)
  xk3 = NewtonMCHomo(as.matrix(trail), zij_ones[,,1:(p+1)], t, h2, n, nn, p+1)
  # xk3 = NewtonMCHomo(as.matrix(trail), zij, t, h2, n, nn, p)
  
  xkkHomo = cbind(xkkHomo, matrix(xk3))
  cat(t, "\n")
}

xkkHomo = xkkHomo[, -1]
######
# Ours
a_b_g_hat <- list()
for (pk in 1:(2 * n + p - 1)) {
  a_b_g_hat[[pk]] <- local({
    current_pk <- pk  # 捕获当前循环的值
    function(s) {
      fit <- smooth.spline(tseq, xkk[current_pk, ])
      return(predict(fit, s)$y)
    }
  })
}

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

es_Ni_list <- list()
for( i in c(1:n)){
  cat(i, "\n")
  trail %>% filter(x == i, time >= 0.1, time < 0.9) -> dd
  if(nrow(dd)==0) next
  dd[order(dd$time), ] -> dd
  dd$N = 1:nrow(dd)
  if(nrow(dd) > 200){
    ind = sample(1:nrow(dd), 200, replace = F) %>% sort()
  }else{
    ind = 1:nrow(dd)
  }
  es_Ni <- c()
  for(t in dd$time[ind]){
    # cat(t, ", ")
    int_exp_i = 0
    for(j in 1:n){
      if(j != i)
        int_exp_i = int_exp_i + int_exp_ij[[paste0(i, "->", j)]](t)
    }
    es_Ni <- c(es_Ni, int_exp_i)
  }
  es_Ni_list[[i]] = rbind(N = dd$N[ind], es_Ni = es_Ni)
}
# plot(dd$N[ind], es_Ni)
save(es_Ni_list, file = "goodness_of_fit_bike2/es_Ni.rdata")


es_Nj_list <- list()
for( j in c(1:n)){
  cat(j, "\n")
  trail %>% filter(y == j, time >= 0.1, time < 0.9) -> dd
  if(nrow(dd)==0) next
  dd[order(dd$time), ] -> dd
  dd$N = 1:nrow(dd)
  if(nrow(dd) > 200){
    ind = sample(1:nrow(dd), 200, replace = F) %>% sort()
  }else{
    ind = 1:nrow(dd)
  }
  es_Nj <- c()
  for(t in dd$time[ind]){
    int_exp_j = 0
    for(i in 1:n){
      if(i != j)
        int_exp_j = int_exp_j + int_exp_ij[[paste0(i, "->", j)]](t)
    }
    es_Nj <- c(es_Nj, int_exp_j)
  }
  es_Nj_list[[j]] = rbind(N = dd$N[ind], es_Nj = es_Nj)
}
save(es_Nj_list, file = "goodness_of_fit_bike2/es_Nj.rdata")




##################
# Homo
homo_hat <- list()
for (pk in 1:(p + 1)) {
  homo_hat[[pk]] <- local({
    current_pk <- pk  # 捕获当前循环的值
    function(s) {
      fit <- smooth.spline(tseq, xkkHomo[current_pk, ])
      return(predict(fit, s)$y)
    }
  })
}


int_exp_ij_homo <- exp_ij_homo <- zij_list <- gs_list <- list()
for (i in 1:n) {
  for (j in 1:n) {
    exp_ij_homo[[paste0(i, "->", j)]] <- local({
      current_i <- i
      current_j <- j
      current_n <- n
      for(pp in 1:(p+1)){
        zij_list[[pp]] <- zij_ones[current_i, current_j, pp]
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
for( i in c(1:n)){
  cat(i, "\n")
  trail %>% filter(x == i, time >= 0.1, time < 0.9) -> dd
  if(nrow(dd)==0) next
  dd[order(dd$time), ] -> dd
  dd$N = 1:nrow(dd)
  if(nrow(dd) > 200){
    ind = sample(1:nrow(dd), 200, replace = F) %>% sort()
  }else{
    ind = 1:nrow(dd)
  }
  es_Ni_homo <- c()
  for(t in dd$time[ind]){
    # cat(t, ",")
    int_exp_i_homo = 0
    for(j in 1:n){
      if(j != i )
        int_exp_i_homo = int_exp_i_homo + int_exp_ij_homo[[paste0(i, "->", j)]](t)
    }
    es_Ni_homo <- c(es_Ni_homo, int_exp_i_homo)
  }
  es_Ni_homo_list[[i]] = rbind(N = dd$N[ind], es_Ni_homo = es_Ni_homo)
}
plot(dd$N[ind], es_Ni_homo)
save(es_Ni_homo_list, file = "goodness_of_fit_bike2/homo_es_Ni.rdata")


es_Nj_homo_list <- list()
for( j in c(1:n)){
  cat(j, "\n")
  trail %>% filter(y == j, time >= 0.1, time < 0.9) -> dd
  if(nrow(dd)==0) next
  dd[order(dd$time), ] -> dd
  dd$N = 1:nrow(dd)
  if(nrow(dd) > 200){
    ind = sample(1:nrow(dd), 200, replace = F) %>% sort()
  }else{
    ind = 1:nrow(dd)
  }
  es_Nj_homo <- c()
  for(t in dd$time[ind]){
    # cat(t, ",")
    int_exp_j_homo = 0
    for(i in 1:n){
      if(i != j )
        int_exp_j_homo = int_exp_j_homo + int_exp_ij_homo[[paste0(i, "->", j)]](t)
    }
    es_Nj_homo <- c(es_Nj_homo, int_exp_j_homo)
  }
  es_Nj_homo_list[[j]] = rbind(N = dd$N[ind], es_Nj_homo = es_Nj_homo)
}
# plot(dd$N[ind], es_Nj_homo)
# for( i in c(1, 26, 51, 76, 100)){
#   plot(es_Ni_homo_list[[i]][1, ],
#        es_Ni_homo_list[[i]][2, ])
# }
save(es_Nj_homo_list, file = "goodness_of_fit_bike2/homo_es_Nj.rdata")


