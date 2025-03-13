rm(list = ls())
library(matrixcalc)
library(Rcpp)
library(RcppEigen)

sourceCpp(file = "MatrixSolverC.cpp")
sourceCpp(file = "NT.cpp")
sourceCpp(file = "NT_homo.cpp")
sourceCpp(file = "SimSet_b.cpp")
sourceCpp(file = "CI.cpp")
fg <- function(t, i) {
  return(sin(2*pi*t)/3)
} 
set.seed(100)
# Test 3 ------------------------------------------------------------------
n = 200
p = 1
t = 0


# b = 0
for(b in c(0, 1/3, 1/2, 1)){
  cat('\n', round(b, 2), ': ')
  xkkk = list()
  xkkk2 = list()
  cont = 1
  nnn = 0
  
  while (cont <= 1000) {
    cat(cont, ', ')
    # Generate
    
    # zij = array(rnorm(n*n*p), c(n,n,p))
    
    zij = array(0, c(n,n,p))
    for (i in 1:(n-1)) {
      if (i < 5)
        zij[i,1:(2*n/6),1] = 1
    }
    
    ones = array(1, c(n, n, 1))
    zij_new = abind::abind(zij, ones, along = 3)
    
    trail_sim = SimSetC(n, b, zij)
    trail_sim = as.data.frame(trail_sim)
    
    rv = order(trail_sim[,3])
    trail_sim = trail_sim[rv,]
    nn = nrow(trail_sim)
    colnames(trail_sim) <- c("s", "r", "t")
    
    Nij = matrix(0, nrow = n, ncol = n)
    for (z in 1:nn) {
      i = trail_sim[z, 1]
      j = trail_sim[z, 2]
      Nij[i, j] = Nij[i, j] + 1
    }
    
    # Newton Method -----------------------------------------------------------
    
    xkk = matrix(rep(0, 2*n + p - 1))
    xkk2 = matrix(rep(0, p+1))
    
    for (t in seq(0.1,0.9,0.05)) {
      
      h1 = 0.1*n^(-0.1)
      h2 = 0.015*n^(-0.2)
      
      xk1 = NewtonMC(as.matrix(trail_sim), zij, t, h1, h2, n, nn, p)
      xk2 = NewtonMCHomo(as.matrix(trail_sim), zij_new, t, h2, n, nn, p+1)
      
      # SD ----------------------------------------------------------------------
      
      xkk = cbind(xkk, matrix(xk1))
      xkk2 = cbind(xkk2, matrix(xk2))
    }
    
    xkk = xkk[, -1]
    xkk2 = xkk2[, -1]
    
    xkkCI = xkk
    
    count = 1
    for (t in seq(0.1,0.9,0.05)) {
      
      xk = xkk[, count]
      
      test = CI(as.matrix(trail_sim), zij, xk, t, h1, h2, n, nn, p) 
      
      xkkCI[, count] = test
      count = count + 1
    }
    
    xkkk[[cont]] = rbind(xkk, xkkCI)
    xkkk2[[cont]] = xkk2
    nnn = nn + nnn
    cont = cont + 1
  }
  
  # xkkk
  # cont = 101
  # xkkkk = array(0, dim = c(400, 17, cont-1))
  # xkkkke = array(0, dim = c(400, 17, cont-1))
  # xkkkksd = array(0, dim = c(400, 17, cont-1))
  
  xkkkk = array(0, dim = c(nrow(xkk), ncol(xkk), cont - 1))
  xkkkke = array(0, dim = c(nrow(xkk), 17, cont - 1))
  xkkkksd = array(0, dim = c(nrow(xkk), 17, cont - 1))
  xkkkkhomo = array(0, dim = c(p+1, 17, cont-1))
  
  for (i in 1:(cont-1)) {
    xkkkke[,,i] = xkkk[[i]][1:(2*n-1+p),]
    xkkkksd[,,i] = xkkk[[i]][(2*n+p):(4*n-2+2*p),]
    xkkkkhomo[,,i] = xkkk2[[i]][1:(p+1),]
  }
  xkkkk = xkkkke
  
  # Plot
  
  library(ggplot2)
  
  # confidence band ----------------------------------------------------------
  
  excl = c()
  for (i in 1:(cont-1)) {
    if (is.nan(sum(xkkkke[,,i]))) {
      excl = c(excl, i)
    }
  }
  if (!is.null(excl)) {
    xkkkke = xkkkke[,,-excl]
  }
  
  xkk_sum = apply(xkkkke, c(1, 2), mean, na.rm = TRUE)
  xkk_sum_sd = apply(xkkkksd, c(1, 2), mean, na.rm = TRUE)
  xl = apply(xkkkke, c(1, 2), quantile, probs = 0.025, na.rm = TRUE)
  xu = apply(xkkkke, c(1, 2), quantile, probs = 0.975, na.rm = TRUE)
  
  xkk_sum_homo = apply(xkkkkhomo, c(1, 2), mean)
  xl_homo = apply(xkkkkhomo, c(1, 2), quantile, probs = 0.05)
  xu_homo = apply(xkkkkhomo, c(1, 2), quantile, probs = 0.95)
  
  pk = 2*n-1+p
  p1 = data.frame(t = rep(seq(0.1,0.9,0.05),2), 
                  y = c(xkk_sum[pk,], xkk_sum_homo[1,]),
                  yl = c(xl[pk,], xl_homo[1,]),
                  yu = c(xu[pk,], xu_homo[1,]),
                  het = c(rep("y", 17), rep("n", 17)))
  
  ggplot(p1, aes(x = t, y = y, group = het, color = het, fill = het)) +
    geom_line(size = 0.7) +
    stat_function(fun = fg, args = list(i = pk), color = "black", size = 0.7) + 
    geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0, linetype = "dashed", size = 0.7) + 
    scale_color_manual(values=c("#619CFF", "red")) +
    xlab(expression(italic("t"))) +
    ylab(expression(hat(italic(gamma)))) +
    theme_bw() +
    theme(legend.position = "none") -> pp
  pdf(file = paste0("simu_fig3/b=",round(b, 2),".pdf"))
  plot(pp)
  dev.off()
}


