rm(list = ls())
library(parallel)
library(doParallel)
library(foreach)
library(matrixcalc)
library(Rcpp)
library(RcppEigen)
library(MASS)
library(ggplot2)

# Sys.setenv("PKG_CPPFLAGS" = "-march=native")

sourceCpp(file = "NT.cpp", verbose = TRUE, rebuild = TRUE)
sourceCpp(file = "SimSet.cpp", verbose = TRUE, rebuild = TRUE)
sourceCpp(file = "CI.cpp", verbose = TRUE, rebuild = TRUE)
genzij1 <- function(n, rho) {
  # Create the covariance matrix
  Sigma <- matrix(0, n, n)
  diag(Sigma) <- 1  # Identity matrix for variances
  Sigma[lower.tri(Sigma)] <- rho  # Set off-diagonal elements
  Sigma[upper.tri(Sigma)] <- rho  # Symmetric matrix
  
  # Generate random vector from multivariate normal distribution
  random_vector <- mvrnorm(1, mu = rep(0, n), Sigma = Sigma)
  
  return(random_vector)
}
genzij2<- function(n, rho) {
  # Create the covariance matrix
  Sigma <- matrix(0, n, n)
  for (s in 1:n) {
    for (t in 1:n) {
      if (s == t) {
        Sigma[s, t] <- 3  # Variance
      } else if (abs(s - t) <= 3) {
        Sigma[s, t] <- rho  # Covariance for |s-t| <= 3
      }
    }
  }
  # Generate random vector from multivariate normal distribution
  random_vector <- mvrnorm(1, mu = rep(0, n), Sigma = Sigma)
  return(random_vector)
}

n_set <- c(100)
rho_set <- c(0.2, 0.4, 0.6, 0.8)

Total_cont = 1000 # total repeats
p = 2 # covariates dimention


for(n in n_set){
  for(rho in rho_set){
    
    cat("n =", n, "; rho =", rho)
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
    starttime <- Sys.time()
    
    folder_name <- paste0("simu_z_dependent/re_n=", n, "rho=", rho)
    if (!dir.exists(folder_name)) {
      dir.create(folder_name)
      # cat("文件夹已创建：", folder_name, "\n")
    } else {
      # cat("文件夹已存在：", folder_name, "\n")
    }
    
    # 4 groups
    xkkk = list()
    cont = 1
    nnn = 0
    cat("; cont: ")
    while (cont <= Total_cont) {
      
      # Generate
      zij = array(0, c(n,n,p));
      for(i in 1:n){
        zij[i,,1]=(genzij1(n,rho)>0);
        zij[i,,2]= genzij2(n,rho);
      }
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
      
      # Newton Method -----------------------------------------------------------
      
      xkk = matrix(rep(0, 2*n + p - 1)) # for all parameters
      # xkk2 = matrix(rep(0, p+1))
      cat(cont, ", ")
      for (t in seq(0.1,0.9,0.05)) {
        # cat(t, ", ")
        
        h1 = 0.1*n^(-0.1)
        h2 = 0.015*n^(-0.2)
        
        xk1 = NewtonMC(as.matrix(trail_sim), zij, t, h1, h2, n, nn, p)
        # xk2 = NewtonMCHomo(as.matrix(trail_sim), zij_new, t, h2, n, nn, p+1)
        
        xkk = cbind(xkk, matrix(xk1))
        # xkk2 = cbind(xkk2, matrix(xk2))
      }
      
      xkk = xkk[, -1]
      # xkk2 = xkk2[, -1]
      
      xkkCI = xkk
      
      count = 1
      for (t in seq(0.1,0.9,0.05)) {
        
        xk = xkk[, count]
        
        test = CI(as.matrix(trail_sim), zij, xk, t, h1, h2, n, nn, p) 
        
        
        
        xkkCI[, count] = test
        count = count + 1
      }
      
      
      xkkk[[cont]] = rbind(xkk, xkkCI)
      # xkkk2[[cont]] = xkk2
      nnn = nn + nnn
      cont = cont + 1
    }
    cont = Total_cont+1
    
    xkkkk = array(0, dim = c(nrow(xkk), ncol(xkk), cont - 1))
    xkkkke = array(0, dim = c(nrow(xkk), 17, cont - 1))
    xkkkksd = array(0, dim = c(nrow(xkk), 17, cont - 1))
    
    for (i in 1:(Total_cont)) {
      xkkkke[,,i] = xkkk[[i]][1:(2*n-1+p),]
      xkkkksd[,,i] = xkkk[[i]][(2*n+p):(4*n-2+2*p),]
    }
    
    
    
    # COVERAGE PROBABILITY, TABEL 2
    library(latex2exp)
    e1 = 0
    pk = c(1, n/2, n+1, n+n/2, 2*n, 2*n+1)
    tk = c(7, 11, 15)
    for (i in 1:(Total_cont)){
      temp = abs(xkkkke[,,i] - X_true)/1.96
      e1 = e1 + (temp[pk, tk] > xkkkksd[pk, tk, i])
    }
    e1 = e1 / (Total_cont)
    1-e1 -> CI_covarages
    
    xkk_sum_sd = apply(xkkkksd, c(1, 2), mean, na.rm = TRUE)
    round(xkk_sum_sd[pk, tk]*1.96*2, 2) -> CI_length
    
    cbind(CI_covarages, CI_length) -> CIs
    colnames(CIs) <- c("t=0.4", "t=0.6", "t=0.8", "t=0.4", "t=0.6", "t=0.8")
    rownames(CIs) <- c("alpha_1", "alpha_n/2+1", "beta_1", "beta_n/2+1", "gamma1", "gamma2")
    CIs
    save(CIs, file = paste0("simu_z_dependent/re_n=", n, "rho=", rho, "/CIs.rdata"))
    
    
    
    # MISE for Table 1--------------------------------------------------------------------
    pk = 1
    fMISE <- function(x, pk, i) {
      fit = smooth.spline(seq(0.1,0.9,0.05), xkkkke[pk,,i])
      return((fs(x, pk) - predict(fit,x)$y)^2)
    }
    
    fMISEr <- function(x, pk, i) {
      fit = smooth.spline(seq(0.1,0.9,0.05), xkkkke[pk + n,,i])
      return((fr(x, pk) - predict(fit,x)$y)^2)
    }
    
    fMISEg <- function(x, i) {
      fit = smooth.spline(seq(0.1,0.9,0.05), xkkkke[2*n,,i])
      return((fg(x, 1) - predict(fit,x)$y)^2)
    }
    
    fMISEg2 <- function(x, i) {
      fit = smooth.spline(seq(0.1,0.9,0.05), xkkkke[2*n+1,,i])
      return((fg(x, 1) - predict(fit,x)$y)^2)
    }
    
    fM = 0
    fM1 = 0
    fM2 = 0
    fM3 = 0
    fM4 = 0
    fM5 = 0
    for (i in 1:Total_cont) {
      for (pk in 1:1) {
        fM = fM + integrate(fMISE, 0, 1, pk = pk, i = i)[[1]]
      }
      
      for (pk in (n/2+1):(n/2+1)) {
        fM1 = fM1 + integrate(fMISE, 0, 1, pk = pk, i = i)[[1]]
      }
      
      for (pk in 1:1) {
        fM2 = fM2 + integrate(fMISEr, 0, 1, pk = pk, i = i)[[1]]
      }
      
      for (pk in (n/2+1):(n/2+1)) {
        fM3 = fM3 + integrate(fMISEr, 0, 1, pk = pk, i = i)[[1]]
      }
      fM4 = fM4 + integrate(fMISEg, 0, 1, i = i)[[1]]
      fM5 = fM5 + integrate(fMISEg2, 0, 1, i = i)[[1]]
    }
    
    
    MISEs = c(fM/Total_cont, fM1/Total_cont, fM2/Total_cont, fM3/Total_cont, fM4/Total_cont, fM5/Total_cont)
    names(MISEs) = c("alpha_1", "alpha_n/2+1", "beta_1", "beta_n/2+1", "gamma1",  "gamma2")
    MISEs
    
    save(MISEs, file = paste0("simu_z_dependent/re_n=", n, "rho=", rho, "/MISES.rdata"))
    
    
    
    # FIG 1 PLOT --------------------------------------------------------------
    excl = c()
    for (i in 1:(Total_cont)) {
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
    
    # PLOT ALPHA
    pk = 1 # pk denote the individual id ranging from 1 to n 
    p1 = data.frame(t = seq(0.1,0.9,0.05), 
                    y = xkk_sum[pk,],
                    yl = xl[pk,],
                    yu = xu[pk,])
    write.csv(p1, file = paste0("simu_z_dependent/re_n=", n, "rho=", rho, "/alpha1.csv"))
    pdf(paste0("simu_z_dependent/re_n=", n, "rho=", rho, "/alpha1.pdf"))
    ggplot(p1, aes(x = t, y = y)) +
      geom_line(color = "red") +
      stat_function(fun = fs, args = list(i = pk), color = "black") + 
      geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0, color = "red", linetype = "dashed") + 
      xlab(expression(italic("t"))) +
      ylab(expression(hat(italic(alpha)))) +
      theme_bw() -> p2
    plot(p2)
    dev.off()
    # PLOT BETA
    pk = 1
    p2 = data.frame(t = seq(0.1,0.9,0.05), 
                    y = xkk_sum[pk+n,],
                    yl = xl[pk+n,],
                    yu = xu[pk+n,])
    write.csv(p1, file = paste0("simu_z_dependent/re_n=", n, "rho=", rho, "/beta1.csv"))
    pdf(paste0("simu_z_dependent/re_n=", n, "rho=", rho, "/beta1.pdf"))
    ggplot(p2, aes(x = t, y = y)) +
      geom_line(color = "red") +
      stat_function(fun = fr, args = list(i = pk), color = "black") + 
      geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0, color = "red", linetype = "dashed") + 
      xlab(expression(italic("t"))) +
      ylab(expression(hat(italic(beta)))) +
      theme_bw() -> p2
    plot(p2)
    dev.off()
    
    # PLOT GAMMA1
    pk = 2*n
    p1 = data.frame(t = rep(seq(0.1,0.9,0.05),1), 
                    y = c(xkk_sum[pk,]),
                    yl = c(xl[pk,]),
                    yu = c(xu[pk,]))
    write.csv(p1, file = paste0("simu_z_dependent/re_n=", n, "rho=", rho, "/gamma1.csv"))
    pdf(paste0("simu_z_dependent/re_n=", n, "rho=", rho, "/gamma1.pdf"))
    ggplot(p1, aes(x = t, y = y)) +
      geom_line(color = "red") +
      stat_function(fun = fg, args = list(i = pk), color = "black") + 
      geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0, color = "red", linetype = "dashed") + 
      xlab(expression(italic("t"))) +
      ylab(expression(hat(italic(gamma)))) +
      theme_bw() -> p2
    plot(p2)
    dev.off()
    
    # PLOT GAMMA2
    pk = 2*n+1
    p1 = data.frame(t = rep(seq(0.1,0.9,0.05),1), 
                    y = c(xkk_sum[pk,]),
                    yl = c(xl[pk,]),
                    yu = c(xu[pk,]))
    write.csv(p1, file = paste0("simu_z_dependent/re_n=", n, "rho=", rho, "/gamma2.csv"))
    
    pdf(paste0("simu_z_dependent/re_n=", n, "rho=", rho, "/gamma2.pdf"))
    ggplot(p1, aes(x = t, y = y)) +
      geom_line(color = "red") +
      stat_function(fun = fg, args = list(i = pk), color = "black") + 
      geom_ribbon(aes(ymin = yl, ymax = yu), alpha = 0, color = "red", linetype = "dashed") + 
      xlab(expression(italic("t"))) +
      ylab(expression(hat(italic(gamma)))) +
      theme_bw() -> p2
    plot(p2)
    dev.off()
    
    endtime <- Sys.time()
    cat(". time:", endtime - starttime, "\n")
  }
}





