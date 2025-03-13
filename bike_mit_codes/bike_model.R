library(RcppEigen)
library(ggplot2)
library(Rcpp)
library(matrixcalc)
library(dplyr)
rm(list = ls())
# setwd("C:/Users/micha/Downloads/dynamic_network")
sourceCpp(file = "cpp/NT_test.cpp")
# sourceCpp(file = "cpp/NT_cv.cpp")
sourceCpp(file = "cpp/NT_ini.cpp")
sourceCpp(file = "cpp/CI.cpp")
sourceCpp(file = "cpp/CI_alt.cpp")
sourceCpp(file = "cpp/NT_homo.cpp")
load(file = "bikedata/bike_data.rdata")

# zij <- bike_data$zij[,,c(3,4)]
# trail <- bike_data$trail

zij <- bike_data$zij[1:417,1:417,c(3,4)]
trail <- bike_data$trail %>% filter(x <= 417, y<= 417)


n = dim(zij)[1]
nn = nrow(trail)
p = 2

Nij = matrix(0, nrow = n, ncol = n)
for (z in 1:nn) {
  i = trail[z, 1] %>% as.numeric()
  j = trail[z, 2] %>% as.numeric()
  Nij[i, j] = Nij[i, j] + 1
}
write.csv(Nij, file = "bikeplots/Nij_bike.csv")
colSums(Nij) -> colsum_bike
rowSums(Nij) -> rowsum_bike
cbind(colsum_bike, rowsum_bike) -> dd_bike
colnames(dd_bike) = c("colsum", "rowsum")
write.csv(dd_bike, file = "bikeplots/dd_bike.csv")

# h1 = 0.025*n^(-0.1)#0.01577393
# h2 = 0.05*n^(-0.2)#0.01990536
# 
# h1 = 0.06*n^(-0.1)#0.01577393
# h2 = 0.12*n^(-0.2)#0.01990536
# 
h1 = 0.08*n^(-0.1)#0.04375991
h2 = 0.16*n^(-0.2)#0.04787323

xkk = matrix(rep(0, 2*n+p-1))

cat(h1, h2, "\n")

xk1 = NewtonMC_ini(as.matrix(trail), zij, rep(0, 2*n + p - 1), 0.1, h1, h2, n, nn, p)
# xk1 = NewtonMC_ini(as.matrix(trail), zij, c(xk1, 0), 0.1, h1, h2, n, nn, 16)

# sum(trail$time <0.1)/nn

for (t in seq(0.1, 0.9, 0.01)) {
  xk2 = NewtonMC_test(as.matrix(trail), zij, rep(0, 2*n + p - 1), t, h1, h2, n, nn, p)
  xkk = cbind(xkk, matrix(xk2))
  cat(t, "\n")
}


xkk = xkk[, -1]

# Plot Summary ------------------------------------------------------------


xkkCI = xkk
count = 1

for (t in seq(0.1, 0.9, 0.01)) {
  
  xk = xkk[, count]
  
  test1 = CI(as.matrix(trail), zij, xk, t, h1, h2, n, nn, p) 
  test =  sqrt(diag((test1$S*1000) %*% (test1$W) %*% (test1$S*1000)/1000000/n^2))
  
  xkkCI[, count] = c(test, test1$xkCI[(2*n):(2*n + p - 1)])
  count = count + 1
}


# # test
# 
# xk = xkk[, 40]
# test = CI(as.matrix(trail), zij, xk, seq(0.1,0.9,0.01)[40], h1, h2, n, nn, p)
# # sqrt(diag(test$O1))
# 
# diag(test$S)
# test$S[1:10,1:10]
# 
# 1/test$vn
# 1/diag(test$V)
# 
# TS = test$S/test$S[1]
# O1 = TS %*% test$W %*% TS
# 
# O1 = (test$S*1000) %*% (test$W) %*% (test$S*1000)/1000000/n^2
# O1 = ((test$S) %*% (test$W) %*% (test$S))/n^2
# 
# #O2 = strassen(strassen(test$S, test$W), test$S)
# isSymmetric(test$S)
# isSymmetric(test$W)
# isSymmetric(O1)
# diag(O1)
# diag(test1$O1)


# Plots -------------------------------------------------------------------

# xkk = xkk[,1:17]
# xkkCI = xkkCI[,1:17]
# seq(0.1,0.9,0.05)[1:17]



ave1 = colSums(xkk[1:n,])/n
ave2 = colSums(xkk[(n+1):(2*n-1),])/(n-1)

pk = 3

ftrans = function(t) {
  return(as.Date("20180501", format = "%Y%m%d") + (as.Date("20180531", format = "%Y%m%d") - as.Date("20180501", format = "%Y%m%d")) * t)
}

# ALPHA
p1 = data.frame(t = ftrans(seq(0.1,0.9,0.01)), 
                y = xkk[pk,] - ave1,
                yl = xkk[pk,] - ave1 - 1.96*xkkCI[pk,],
                yu = xkk[pk,] - ave1 + 1.96*xkkCI[pk,])

# BETA
p2 = data.frame(t = ftrans(seq(0.1,0.9,0.01)), 
                y = xkk[pk+n,] - ave2,
                yl = xkk[pk+n,] - ave2 - 1.96*xkkCI[pk,],
                yu = xkk[pk+n,] - ave2 + 1.96*xkkCI[pk,])

# 1a
ggplot(p1, aes(x = t, y = y)) +
  geom_line(color = "red", size = 0.75) +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0, color = "red", linetype = "dashed", size = 0.75) +
  scale_x_date(date_labels = "%b %d") +
  xlab(expression(italic("t"))) +
  ylab(expression(widehat(italic(alpha))[i]^"'"~(t))) +
  ylim(-5, 5)+
  theme(panel.background = element_rect(fill = "white")) + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

# 1b
ggplot(p2, aes(x = t, y = y)) +
  geom_line(color = "red", size = 0.75) +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0, color = "red", linetype = "dashed", size = 0.75) +
  scale_x_date(date_labels = "%b %d") +
  xlab(expression(italic("t"))) +
  ylab(expression(widehat(italic(beta))[i]^"'"~(t))) +
  ylim(-5, 5)+
  theme(panel.background = element_rect(fill = "white")) + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))




# xkkHomo = matrix(rep(0, p))
xkkHomo = matrix(rep(0, p+1))
ones = array(1, c(n, n, 1))
zij_ones = abind::abind(zij, ones, along = 3)


for (t in seq(0.1,0.9,0.01)) {
  # xk2 = NewtonMC_test(as.matrix(trail), zij, xk1, t, h1, h2, n, nn, 16)
  # xk3 = NewtonMCHomo(as.matrix(trail), zij_ones[,,5:6], t, h2, n, nn, p+1)
  # xk3 = NewtonMCHomo(as.matrix(trail), zij_ones[,,c(1,6)], t, h2, n, nn, 5)
  xk3 = NewtonMCHomo(as.matrix(trail), zij_ones[,,1:(p+1)], t, h2, n, nn, p+1)
  # xk3 = NewtonMCHomo(as.matrix(trail), zij, t, h2, n, nn, p)
  
  xkkHomo = cbind(xkkHomo, matrix(xk3))
  cat(t, "\n")
}

rerun = which(is.nan(xkkHomo[1,]))

for (i in rerun) {
  t = seq(0.1,0.9,0.01)[i-1]
  # xk2 = NewtonMC_test(as.matrix(trail), zij, xk1, t, h1, h2, n, nn, 16)
  # xk3 = NewtonMCHomo(as.matrix(trail), zij_ones[,,5:6], t, h2, n, nn, p+1)
  # xk3 = NewtonMCHomo(as.matrix(trail), zij_ones[,,c(1,6)], t, h2, n, nn, 5)
  xk3 = NewtonMCHomo(as.matrix(trail), zij_ones[,,1:6], t, h2, n, nn, p+1)
  # xk3 = NewtonMCHomo(as.matrix(trail), zij, t, h2, n, nn, p)
  
  xkkHomo[,i] = xk3
  cat(t, "\n")
}

# xkk = xkk[, -1]
xkkHomo = xkkHomo[, -1]

save(xkk, file = paste0("bikedata/xkk", h1,"_1.rdata"))
save(xkkCI, file = paste0("bikedata/xkkCI", h1,"_1.rdata"))
save(xkkHomo, file = paste0("bikedata/xkkHomo", h1,"_1.rdata"))


# plot(xkk[,33], xkk[,34])
# plot(xkkCI[1:300,50], xkkCI[1:300,51])
# plot(xkkCI[,50], xkkCI[,51])

# xk1 = xkk[, 33]
# t = t1[33]
# test1 = CI(as.matrix(trail), zij, xk1, t, h1, h2, n, nn, p) 

# xk2 = xkk[, 34]
# t = t1[34]
# test2 = CI(as.matrix(trail), zij, xk2, t, h1, h2, n, nn, p) 


# plot(test1, test2)
# plot(xk1, xk2)

pp = 2 # change according to pk (1,2,3,4,5)
pk = 2*n-1+pp

t2 = seq(0.1,0.9,0.01)

pp=1
p1 = data.frame(t = seq(0.1,0.9,0.01)%>%ftrans(), 
                y = c(xkk[2*n - 1 + pp,], xkkHomo[pp,]),
                yl = c(xkk[2*n - 1 + pp,] - 1.96*xkkCI[2*n - 1 + pp,], xkkHomo[pp,]),
                yu = c(xkk[2*n - 1 + pp,] + 1.96*xkkCI[2*n - 1 + pp,], xkkHomo[pp,]),
                het = c(rep("y", length(t2)), rep("n", length(t2))))

pp=2
p2 = data.frame(t = seq(0.1,0.9,0.01)%>%ftrans(), 
                y = c(xkk[2*n - 1 + pp,], xkkHomo[pp,]),
                yl = c(xkk[2*n - 1 + pp,] - 1.96*xkkCI[2*n - 1 + pp,], xkkHomo[pp,]),
                yu = c(xkk[2*n - 1 + pp,] + 1.96*xkkCI[2*n - 1 + pp,], xkkHomo[pp,]),
                het = c(rep("y", length(t2)), rep("n", length(t2))))
# sf

ggplot(p1, aes(x = t, y = y, group = het, colour = het, fill = het)) +
  geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0, linetype = "dashed", size = 0.7) + 
  scale_color_manual(values=c("#619CFF", "red")) +
  scale_x_date(date_labels = "%b %d", breaks = "6 week") +
  xlab(expression(italic("t"))) +
  ylab(expression(widehat(italic(gamma))[1](t))) +
  theme(panel.background = element_rect(fill = "white")) + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(legend.position = "none")

# s0

ggplot(p2, aes(x = t, y = y, group = het, colour = het, fill = het)) +
  geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0, linetype = "dashed", size = 0.7) + 
  scale_color_manual(values=c("#619CFF", "red")) +
  xlab(expression(italic("t"))) +
  ylab(expression(widehat(italic(gamma))[2](t))) +
  theme(panel.background = element_rect(fill = "white")) + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(legend.position = "none")

###################
# plot bike results
library(ggplot2)
library(dplyr)

load(file = "bikedata/xkk0.0437599061077816_1.rdata")
load(file = "bikedata/xkkCI0.0437599061077816_1.rdata")
load(file = "bikedata/xkkHomo0.0437599061077816_1.rdata")
ftrans = function(t) {
  return(as.Date("20180501", format = "%Y%m%d") + (as.Date("20180531", format = "%Y%m%d") - as.Date("20180501", format = "%Y%m%d")) * t)
}
n = 417
ave1 = colSums(xkk[1:n,])/n
ave2 = colSums(xkk[(n+1):(2*n-1),])/(n-1)
for(pk in c(1, 101, 201, 301, 401)){
  # ALPHA
  p1 = data.frame(t1 = seq(0.1,0.9,0.01) %>% ftrans(), 
                  t2 = seq(0.1, 0.9, 0.01), 
                  y = xkk[pk,] - ave1,
                  yl = xkk[pk,] - ave1 - 1.96*xkkCI[pk,],
                  yu = xkk[pk,] - ave1 + 1.96*xkkCI[pk,])
  
  # BETA
  p2 = data.frame(t1 = seq(0.1,0.9,0.01) %>% ftrans(), 
                  t2 = seq(0.1, 0.9, 0.01),  
                  y = xkk[pk+n,] - ave2,
                  yl = xkk[pk+n,] - ave2 - 1.96*xkkCI[pk,],
                  yu = xkk[pk+n,] - ave2 + 1.96*xkkCI[pk,])
  
  # 1a
  ggplot(p1, aes(x = t2, y = y)) +
    geom_line(color = "red", size = 0.75) +
    geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0, color = "red", linetype = "dashed", size = 0.75) +  # 添加fill颜色映射
    scale_x_continuous(
      breaks = seq(0.1, 0.9, by = 0.2),  # 调整间隔为每0.2显示一个标签
      labels = function(t) {
        format(ftrans(t), "%b %d")       # 将数值转换为日期格式
      }
    ) +
    xlab(expression(italic("t"))) +
    ylab(expression(widehat(italic(alpha))[i]^"'"~(t))) +
    ylim(-5, 5)+
    theme(panel.background = element_rect(fill = "white")) + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) -> f1
  pdf(file = paste0("bikeplots/alpha_", pk, ".pdf"))
  plot(f1)
  dev.off()
  # 1b
  ggplot(p2, aes(x = t2, y = y)) +
    geom_line(color = "red", size = 0.75) +
    geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0, color = "red", linetype = "dashed", size = 0.75) +# 添加fill颜色映射
    scale_x_continuous(
      breaks = seq(0.1, 0.9, by = 0.2),  # 调整间隔为每0.2显示一个标签
      labels = function(t) {
        format(ftrans(t), "%b %d")       # 将数值转换为日期格式
      }
    ) +
    xlab(expression(italic("t"))) +
    ylab(expression(widehat(italic(beta))[i]^"'"~(t))) +
    ylim(-5, 5)+
    theme(panel.background = element_rect(fill = "white")) + 
    theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) -> f2
  pdf(file = paste0("bikeplots/beta_", pk, ".pdf"))
  plot(f2)
  dev.off()
}
