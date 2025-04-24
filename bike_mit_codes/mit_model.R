rm(list = ls())
library(RcppEigen)
library(Rcpp)
library(matrixcalc)

# setwd("C:/Users/micha/Downloads/dynamic_network")
sourceCpp(file = "cpp/NT_test.cpp")
sourceCpp(file = "cpp/NT_ini.cpp")
sourceCpp(file = "cpp/CI.cpp")
sourceCpp(file = "cpp/CI_alt.cpp")
sourceCpp(file = "cpp/NT_homo.cpp")

load(file = "proxyData.rda")

n = proxyData$n
nn = nrow(proxyData$trail)
p = 5
trail = proxyData$trail

Nij = matrix(0, nrow = n, ncol = n)
for (z in 1:nn) {
  i = trail[z, 1]
  j = trail[z, 2]
  Nij[i, j] = Nij[i, j] + 1
}
write.csv(Nij, file = "mit_plots/Nij_mit.csv")
colSums(Nij) -> colsum_MIT
rowSums(Nij) -> rowsum_MIT
cbind(colsum_MIT, rowsum_MIT) -> dd_MIT
colnames(dd_MIT) = c("colsum", "rowsum")
write.csv(dd_MIT, file = "mit_plots/dd_mit.csv")



h1 = 0.1*n^(-0.1)
h2 = 0.015*n^(-0.2)
xkk = matrix(rep(0, 2*n+p-1))

h1 = 0.05
h2 = 0.05

zij = proxyData$zij[,,c(6:10),1]

# xk1 = NewtonMC_ini(as.matrix(trail), array(zij[,,c(1)], c(n, n, 1)), rep(0, 2*n + 1 - 1), 0.1, h1, h2, n, nn, 1)
xk1 = NewtonMC_ini(as.matrix(trail), zij, rep(0, 2*n + p - 1), 0.1, h1, h2, n, nn, p)
# xk1 = NewtonMC_ini(as.matrix(trail), zij, c(xk1, 0), 0.1, h1, h2, n, nn, 16)

sum(trail$time < 0.2)/nn

for (t in seq(0.2,0.9,0.01)) {
  xk2 = NewtonMC_test(as.matrix(trail), zij, xk1, t, h1, h2, n, nn, p)
  xkk = cbind(xkk, matrix(xk2))
  cat(t, "\n")
}

xkk = xkk[, -1]



# Plot Summary ------------------------------------------------------------


xkkCI = xkk
count = 1

for (t in seq(0.2,0.9,0.01)) {
  
  xk = xkk[, count]
  
  test1 = CI(as.matrix(trail), zij, xk, t, h1, h2, n, nn, p) 
  test =  sqrt(diag((test1$S*1000) %*% (test1$W) %*% (test1$S*1000)/1000000/n^2))
  
  xkkCI[, count] = c(test, test1$xkCI[126:130])
  count = count + 1
}


# test

xk = xkk[, 40]
test = CI(as.matrix(trail), zij, xk, 0.59, h1, h2, n, nn, p) 
sqrt(diag(test$O1))

diag(test$S)
test$S[1:10,1:10]

1/test$vn
1/diag(test$V)

TS = test$S/test$S[1]
O1 = TS %*% test$W %*% TS

O1 = (test$S*1000) %*% (test$W) %*% (test$S*1000)/1000000/n^2
O1 = ((test$S) %*% (test$W) %*% (test$S))/n^2

#O2 = strassen(strassen(test$S, test$W), test$S)
isSymmetric(test$S)
isSymmetric(test$W)
isSymmetric(O1)
diag(O1)
diag(test1$O1)


# Plots -------------------------------------------------------------------

xkk = xkk[,5:71]
xkkCI = xkkCI[,5:71]
seq(0.2,0.9,0.01)[5:71]
# save(xkk, file = "mit_plots/xkk.rdata")
# save(xkkCI, file = "mit_plots/xkkCI.rdata")

library(ggplot2)

ave1 = colSums(xkk[1:n,])/n
ave2 = colSums(xkk[(n+1):(2*n-1),])/(n-1)


#  August 24, 1998 -> June 21, 2002
ftrans = function(t) {
  return(as.Date("20080908", format = "%Y%m%d") + (as.Date("20090625", format = "%Y%m%d") - as.Date("20080908", format = "%Y%m%d")) * t)
}

# ALPHA
pk=1
p1 = data.frame(t = ftrans(seq(0.24,0.9,0.01)), 
                y = xkk[pk,] - ave1,
                yl = xkk[pk,] - ave1 - 1.96*xkkCI[pk,],
                yu = xkk[pk,] - ave1 + 1.96*xkkCI[pk,])

# BETA
p2 = data.frame(t = ftrans(seq(0.24,0.9,0.01)), 
                y = xkk[pk+n,] - ave2,
                yl = xkk[pk+n,] - ave2 - 1.96*xkkCI[pk,],
                yu = xkk[pk+n,] - ave2 + 1.96*xkkCI[pk,])

# 1a
ggplot(p1, aes(x = t, y = y)) +
  geom_line(color = "red", size = 0.75) +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0, color = "red", linetype = "dashed", size = 0.75) +
  scale_x_date(date_labels = "%b %d", breaks = "6 week") +
  scale_y_continuous(breaks=seq(-2, 6, 2), limits = c(-2, 6)) +
  xlab(expression(italic("t"))) +
  ylab(expression(widehat(italic(alpha))[i]^"'"~(t))) +
  theme(panel.background = element_rect(fill = "white")) + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

# 1b
ggplot(p2, aes(x = t, y = y)) +
  geom_line(color = "red", size = 0.75) +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0, color = "red", linetype = "dashed", size = 0.75) +
  scale_x_date(date_labels = "%b %d", breaks = "6 week") +
  scale_y_continuous(breaks=seq(-0, 5, 1), limits = c(-0, 5)) +
  xlab(expression(italic("t"))) +
  ylab(expression(widehat(italic(beta))[i]^"'"~(t))) +
  theme(panel.background = element_rect(fill = "white")) + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))



# ALPHA
pk=28
p1 = data.frame(t = ftrans(seq(0.24,0.9,0.01)), 
                y = xkk[pk,] - ave1,
                yl = xkk[pk,] - ave1 - 1.96*xkkCI[pk,],
                yu = xkk[pk,] - ave1 + 1.96*xkkCI[pk,])

# BETA
p2 = data.frame(t = ftrans(seq(0.24,0.9,0.01)), 
                y = xkk[pk+n,] - ave2,
                yl = xkk[pk+n,] - ave2 - 1.96*xkkCI[pk,],
                yu = xkk[pk+n,] - ave2 + 1.96*xkkCI[pk,])

# 28 a
ggplot(p1, aes(x = t, y = y)) +
  geom_line(color = "red", size = 0.75) +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0, color = "red", linetype = "dashed", size = 0.75) +
  scale_x_date(date_labels = "%b %d", breaks = "6 week") +
  scale_y_continuous(breaks=seq(-8, 6, 2), limits = c(-8, 6)) +
  xlab(expression(italic("t"))) +
  ylab(expression(widehat(italic(alpha))[i]^"'"~(t))) +
  theme(panel.background = element_rect(fill = "white")) + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))

# 28 b
ggplot(p2, aes(x = t, y = y)) +
  geom_line(color = "red", size = 0.75) +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0, color = "red", linetype = "dashed", size = 0.75) +
  scale_x_date(date_labels = "%b %d", breaks = "6 week") +
  scale_y_continuous(breaks=seq(-2, 2, 1), limits = c(-2, 2)) +
  xlab(expression(italic("t"))) +
  ylab(expression(widehat(italic(beta))[i]^"'"~(t))) +
  theme(panel.background = element_rect(fill = "white")) + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1))


# GAMMA



# xkkHomo = matrix(rep(0, p))
xkkHomo = matrix(rep(0, p+1))
ones = array(1, c(n, n, 1))
zij_ones = abind::abind(zij, ones, along = 3)


for (t in seq(0.24,0.9,0.01)) {
  # xk2 = NewtonMC_test(as.matrix(trail), zij, xk1, t, h1, h2, n, nn, 16)
  # xk3 = NewtonMCHomo(as.matrix(trail), zij_ones[,,5:6], t, h2, n, nn, p+1)
  # xk3 = NewtonMCHomo(as.matrix(trail), zij_ones[,,c(1,6)], t, h2, n, nn, 5)
  xk3 = NewtonMCHomo(as.matrix(trail), zij_ones[,,1:6], t, h2, n, nn, p+1)
  # xk3 = NewtonMCHomo(as.matrix(trail), zij, t, h2, n, nn, p)
  
  xkkHomo = cbind(xkkHomo, matrix(xk3))
  cat(t, "\n")
}

rerun = which(is.nan(xkkHomo[1,]))

for (i in rerun) {
  t = seq(0.24,0.9,0.01)[i-1]
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
save(xkkHomo, file = "mit_plots/xkkHomo.rdata")
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

pk = 125+1 # same floor
pk = 125+2 # 0
pk = 125+3 # 1
pk = 125+4 # 2
pk = 125+5 # 3


t2 = seq(0.24,0.9,0.01)

p1 = data.frame(t = ftrans(t2), 
                y = xkk[pk,],
                yl = xkk[pk,] - 1.96*xkkCI[pk,],
                yu = xkk[pk,] + 1.96*xkkCI[pk,])

pp = 1 # change according to pk (1,2,3,4,5)
p1 = data.frame(t = ftrans(seq(0.24,0.9,0.01)), 
                y = c(xkk[2*n - 1 + pp,], xkkHomo[pp,]),
                yl = c(xkk[2*n - 1 + pp,] - 1.96*xkkCI[2*n - 1 + pp,], xkkHomo[pp,]),
                yu = c(xkk[2*n - 1 + pp,] + 1.96*xkkCI[2*n - 1 + pp,], xkkHomo[pp,]),
                het = c(rep("y", 67), rep("n", 67)))
ggplot(p1, aes(x = t, y = y, group = het, colour = het, fill = het)) +
  geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0, linetype = "dashed", size = 0.7) + 
  scale_color_manual(values=c("#619CFF", "red")) +
  scale_x_date(date_labels = "%b %d", breaks = "6 week") +
  xlab(expression(italic("t"))) +
  ylab(expression(widehat(italic(gamma))[1](t))) +
  theme(panel.background = element_rect(fill = "white")) + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(legend.position = "none") -> p1

pp = 2 # change according to pk (1,2,3,4,5)
p1 = data.frame(t = ftrans(seq(0.24,0.9,0.01)), 
                y = c(xkk[2*n - 1 + pp,], xkkHomo[pp,]),
                yl = c(xkk[2*n - 1 + pp,] - 1.96*xkkCI[2*n - 1 + pp,], xkkHomo[pp,]),
                yu = c(xkk[2*n - 1 + pp,] + 1.96*xkkCI[2*n - 1 + pp,], xkkHomo[pp,]),
                het = c(rep("y", 67), rep("n", 67)))
ggplot(p1, aes(x = t, y = y, group = het, colour = het, fill = het)) +
  geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0, linetype = "dashed", size = 0.7) + 
  scale_color_manual(values=c("#619CFF", "red")) +
  scale_x_date(date_labels = "%b %d", breaks = "6 week") +
  xlab(expression(italic("t"))) +
  ylab(expression(widehat(italic(gamma))[2](t))) +
  theme(panel.background = element_rect(fill = "white")) + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(legend.position = "none") -> p2



pp = 3 # change according to pk (1,2,3,4,5)
p1 = data.frame(t = ftrans(seq(0.24,0.9,0.01)), 
                y = c(xkk[2*n - 1 + pp,], xkkHomo[pp,]),
                yl = c(xkk[2*n - 1 + pp,] - 1.96*xkkCI[2*n - 1 + pp,], xkkHomo[pp,]),
                yu = c(xkk[2*n - 1 + pp,] + 1.96*xkkCI[2*n - 1 + pp,], xkkHomo[pp,]),
                het = c(rep("y", 67), rep("n", 67)))
ggplot(p1, aes(x = t, y = y, group = het, colour = het, fill = het)) +
  geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0, linetype = "dashed", size = 0.7) + 
  scale_color_manual(values=c("#619CFF", "red")) +
  scale_x_date(date_labels = "%b %d", breaks = "6 week") +
  xlab(expression(italic("t"))) +
  ylab(expression(widehat(italic(gamma))[3](t))) +
  theme(panel.background = element_rect(fill = "white")) + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(legend.position = "none") -> p3


pp = 4 # change according to pk (1,2,3,4,5)
p1 = data.frame(t = ftrans(seq(0.24,0.9,0.01)), 
                y = c(xkk[2*n - 1 + pp,], xkkHomo[pp,]),
                yl = c(xkk[2*n - 1 + pp,] - 1.96*xkkCI[2*n - 1 + pp,], xkkHomo[pp,]),
                yu = c(xkk[2*n - 1 + pp,] + 1.96*xkkCI[2*n - 1 + pp,], xkkHomo[pp,]),
                het = c(rep("y", 67), rep("n", 67)))
ggplot(p1, aes(x = t, y = y, group = het, colour = het, fill = het)) +
  geom_line() +
  geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0, linetype = "dashed", size = 0.7) + 
  scale_color_manual(values=c("#619CFF", "red")) +
  scale_x_date(date_labels = "%b %d", breaks = "6 week") +
  xlab(expression(italic("t"))) +
  ylab(expression(widehat(italic(gamma))[4](t))) +
  theme(panel.background = element_rect(fill = "white")) + 
  theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
  theme(legend.position = "none") -> p4

pdf(file = "mit_plots/gamma1.pdf")
plot(p1)
dev.off()
pdf(file = "mit_plots/gamma2.pdf")
plot(p2)
dev.off()
pdf(file = "mit_plots/gamma3.pdf")
plot(p3)
dev.off()
pdf(file = "mit_plots/gamma4.pdf")
plot(p4)
dev.off()
