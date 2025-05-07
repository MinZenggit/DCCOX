library(RcppEigen)
library(ggplot2)
library(Rcpp)
library(matrixcalc)
library(dplyr)
rm(list = ls())
setwd("~/DCCOX")
sourceCpp(file = "cpp/NT_test.cpp")
# sourceCpp(file = "cpp/NT_cv.cpp")
sourceCpp(file = "cpp/NT_ini.cpp")
sourceCpp(file = "cpp/CI.cpp")
# sourceCpp(file = "cpp/CI_alt.cpp")
sourceCpp(file = "cpp/NT_homo.cpp")
load(file = "bikedata/bike_data3.rdata")

# 
# zij <- bike_data$zij[1:50,1:50,]
# trail <- bike_data$trail %>% filter(x <= 50, y<= 50)

nodes <- 1:542
zij <- bike_data$zij[nodes,nodes,c(1:4)]
zij[is.na(zij)] = 0
trail <- bike_data$trail %>% filter(x %in% nodes, y %in% nodes)
n = dim(zij)[1]
nn = nrow(trail)
p = dim(zij)[3]

Nij = matrix(0, nrow = n, ncol = n)
for (z in 1:nn) {
  i = trail[z, 1]
  j = trail[z, 2]
  Nij[i, j] = Nij[i, j] + 1
}
write.csv(Nij, file = "bikeplots3/Nij_mit.csv")
colSums(Nij) -> colsum_MIT
rowSums(Nij) -> rowsum_MIT
cbind(colsum_MIT, rowsum_MIT) -> dd_MIT
colnames(dd_MIT) = c("colsum", "rowsum")
write.csv(dd_MIT, file = "bikeplots3/dd_mit.csv")

h1 = 0.15*n^(-0.1)#0.04375991
h2 = 0.4*n^(-0.2)#0.04787323
tseq = seq(0.05,0.95,0.01)
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

# Plot Summary ------------------------------------------------------------


xkkCI = xkk
count = 1

for (t in tseq) {
  
  xk = xkk[, count]
  
  test1 = CI(as.matrix(trail), zij, xk, t, h1, h2, n, nn, p) 
  xkkCI[, count] = test1
  count = count + 1
}


ave1 = colSums(xkk[1:n,])/n
ave2 = colSums(xkk[(n+1):(2*n-1),])/(n-1)

pk = 3

ftrans = function(t) {
  return(as.Date("20180101", format = "%Y%m%d") + (as.Date("20181231", format = "%Y%m%d") - as.Date("20180101", format = "%Y%m%d")) * t)
}

# ALPHA
p1 = data.frame(t = ftrans(tseq), 
                y = xkk[pk,] - ave1,
                yl = xkk[pk,] - ave1 - 1.96*xkkCI[pk,],
                yu = xkk[pk,] - ave1 + 1.96*xkkCI[pk,])

# BETA
p2 = data.frame(t = ftrans(tseq), 
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
  ylim(-10, 10)+
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


for (t in tseq) {
  # xk2 = NewtonMC_test(as.matrix(trail), zij, xk1, t, h1, h2, n, nn, 16)
  # xk3 = NewtonMCHomo(as.matrix(trail), zij_ones[,,5:6], t, h2, n, nn, p+1)
  # xk3 = NewtonMCHomo(as.matrix(trail), zij_ones[,,c(1,6)], t, h2, n, nn, 5)
  xk3 = NewtonMCHomo(as.matrix(trail), zij_ones[,,1:(p+1)], t, h2/10, n, nn, p+1)
  # xk3 = NewtonMCHomo(as.matrix(trail), zij, t, h2, n, nn, p)
  
  xkkHomo = cbind(xkkHomo, matrix(xk3))
  cat(t, "\n")
}

xkkHomo = xkkHomo[, -1]

save(xkk, file = paste0("bikeplots3/xkk", h1,"_3.rdata"))
save(xkkCI, file = paste0("bikeplots3/xkkCI", h1,"_3.rdata"))
save(xkkHomo, file = paste0("bikeplots3/xkkHomo", h1,"_3.rdata"))

# 
# # plot(xkk[,33], xkk[,34])
# # plot(xkkCI[1:300,50], xkkCI[1:300,51])
# # plot(xkkCI[,50], xkkCI[,51])
# 
# # xk1 = xkk[, 33]
# # t = t1[33]
# # test1 = CI(as.matrix(trail), zij, xk1, t, h1, h2, n, nn, p) 
# 
# # xk2 = xkk[, 34]
# # t = t1[34]
# # test2 = CI(as.matrix(trail), zij, xk2, t, h1, h2, n, nn, p) 
# 
# 
# # plot(test1, test2)
# # plot(xk1, xk2)
# 
# pp = 2 # change according to pk (1,2,3,4,5)
# pk = 2*n-1+pp
# 
# t2 = tseq
# 
# pp=1
# p1 = data.frame(t = tseq%>%ftrans(), 
#                 y = c(xkk[2*n - 1 + pp,], xkkHomo[pp,]),
#                 yl = c(xkk[2*n - 1 + pp,] - 1.96*xkkCI[2*n - 1 + pp,], xkkHomo[pp,]),
#                 yu = c(xkk[2*n - 1 + pp,] + 1.96*xkkCI[2*n - 1 + pp,], xkkHomo[pp,]),
#                 het = c(rep("y", length(t2)), rep("n", length(t2))))
# 
# pp=2
# p2 = data.frame(t = tseq%>%ftrans(), 
#                 y = c(xkk[2*n - 1 + pp,], xkkHomo[pp,]),
#                 yl = c(xkk[2*n - 1 + pp,] - 1.96*xkkCI[2*n - 1 + pp,], xkkHomo[pp,]),
#                 yu = c(xkk[2*n - 1 + pp,] + 1.96*xkkCI[2*n - 1 + pp,], xkkHomo[pp,]),
#                 het = c(rep("y", length(t2)), rep("n", length(t2))))
# 
# pp=3
# p3 = data.frame(t = tseq%>%ftrans(), 
#                 y = c(xkk[2*n - 1 + pp,], xkkHomo[pp,]),
#                 yl = c(xkk[2*n - 1 + pp,] - 1.96*xkkCI[2*n - 1 + pp,], xkkHomo[pp,]),
#                 yu = c(xkk[2*n - 1 + pp,] + 1.96*xkkCI[2*n - 1 + pp,], xkkHomo[pp,]),
#                 het = c(rep("y", length(t2)), rep("n", length(t2))))
# 
# pp=4
# p4 = data.frame(t = tseq%>%ftrans(), 
#                 y = c(xkk[2*n - 1 + pp,], xkkHomo[pp,]),
#                 yl = c(xkk[2*n - 1 + pp,] - 1.96*xkkCI[2*n - 1 + pp,], xkkHomo[pp,]),
#                 yu = c(xkk[2*n - 1 + pp,] + 1.96*xkkCI[2*n - 1 + pp,], xkkHomo[pp,]),
#                 het = c(rep("y", length(t2)), rep("n", length(t2))))
# # sf
# 
# ggplot(p1, aes(x = t, y = y, group = het, colour = het, fill = het)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0, linetype = "dashed", size = 0.7) + 
#   scale_color_manual(values=c("#619CFF", "red")) +
#   scale_x_date(date_labels = "%b %d", breaks = "6 week") +
#   xlab(expression(italic("t"))) +
#   ylab(expression(widehat(italic(gamma))[1](t))) +
#   theme(panel.background = element_rect(fill = "white")) + 
#   theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
#   theme(legend.position = "none")
# 
# # s0
# 
# ggplot(p2, aes(x = t, y = y, group = het, colour = het, fill = het)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0, linetype = "dashed", size = 0.7) + 
#   scale_color_manual(values=c("#619CFF", "red")) +
#   xlab(expression(italic("t"))) +
#   ylab(expression(widehat(italic(gamma))[2](t))) +
#   theme(panel.background = element_rect(fill = "white")) + 
#   theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
#   theme(legend.position = "none")
# 
# 
# ggplot(p3, aes(x = t, y = y, group = het, colour = het, fill = het)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0, linetype = "dashed", size = 0.7) + 
#   scale_color_manual(values=c("#619CFF", "red")) +
#   xlab(expression(italic("t"))) +
#   ylab(expression(widehat(italic(gamma))[2](t))) +
#   theme(panel.background = element_rect(fill = "white")) + 
#   theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
#   theme(legend.position = "none")
# 
# 
# ggplot(p4, aes(x = t, y = y, group = het, colour = het, fill = het)) +
#   geom_line() +
#   geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0, linetype = "dashed", size = 0.7) + 
#   scale_color_manual(values=c("#619CFF", "red")) +
#   xlab(expression(italic("t"))) +
#   ylab(expression(widehat(italic(gamma))[2](t))) +
#   theme(panel.background = element_rect(fill = "white")) + 
#   theme(panel.border = element_rect(color = "black", fill = NA, size = 1)) +
#   theme(legend.position = "none")
# 
# ###################
# plot bike results
library(ggplot2)
library(dplyr)

load(file = "bikeplots3/xkk0.0532843989577126_2.rdata")
load(file = "bikeplots3/xkkCI0.0532843989577126_2.rdata")
load(file = "bikeplots3/xkkHomo0.0532843989577126_2.rdata")
xkk <- xkk[, 21:91]
xkkCI <- xkkCI[, 21:91]
xkkHomo <- xkkHomo[, 21:91]
tseq = seq(0.05,0.95,0.01)[21:91]
ftrans = function(t) {
  return(as.Date("20180101", format = "%Y%m%d") + (as.Date("20181231", format = "%Y%m%d") - as.Date("20180101", format = "%Y%m%d")) * t)
}
key_dates = c(0.40, 0.68, 0.86)

n = 542
ave1 = colSums(xkk[1:n,])/n
ave2 = colSums(xkk[(n+1):(2*n-1),])/(n-1)
for(pk in c(1, 101, 201, 301, 401, 501)){
  # ALPHA
  p1 = data.frame(t1 = tseq %>% ftrans(),
                  t2 = tseq,
                  y = xkk[pk,] - ave1,
                  yl = xkk[pk,] - ave1 - 1.96*xkkCI[pk,],
                  yu = xkk[pk,] - ave1 + 1.96*xkkCI[pk,])

  # BETA
  p2 = data.frame(t1 = tseq %>% ftrans(),
                  t2 = tseq,
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
    geom_vline(
      xintercept = key_dates, 
      linetype = "dotted", 
      color = "black", 
      size = 1.5,
      alpha = 0.7
    )+
    xlab(expression(italic("t"))) +
    ylab(expression(widehat(italic(alpha))[i]^"'"~(t))) -> f1
  pdf(file = paste0("bikeplots3/alpha_", pk, ".pdf"))
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
    geom_vline(
      xintercept = key_dates, 
      linetype = "dotted", 
      color = "black", 
      size = 1.5,
      alpha = 0.7
    )+
    xlab(expression(italic("t"))) +
    ylab(expression(widehat(italic(beta))[i]^"'"~(t)))  -> f2
  pdf(file = paste0("bikeplots3/beta_", pk, ".pdf"))
  plot(f2)
  dev.off()
}
