NewtonMCHomo_R <- function(Trail, Zij, t, h2, n, nn, p) {
  # 检查Zij的维度
  if (is.null(attr(Zij, "dim"))) {
    stop("'Zij' does not have 'dim' attribute.")
  }
  d <- dim(Zij)
  if (length(d) != 3) {
    stop("'Zij' must have 3 dimensions.")
  }
  if (d[1] != n || d[2] != n || d[3] != p) {
    stop("Dimensions of 'Zij' do not match n, n, p.")
  }
  
  # 初始化变量
  NijT2 <- matrix(0, nrow = n, ncol = n)
  f1 <- numeric(p)
  xk <- numeric(p)
  conv <- 1.0
  iter <- 0
  max_iter <- 1000  # 防止无限循环
  
  # 计算带宽
  bw3 <- max(t - 10*h2, 0)
  bw4 <- min(t + 10*h2, 1)
  
  # 计算NijT2
  inv_sqrt_2pi <- 1 / sqrt(2*pi)
  for (z in 1:nn) {
    time <- Trail[z, 3]
    if (time > bw3 && time < bw4) {
      i <- Trail[z, 1]
      j <- Trail[z, 2]
      diff <- (time - t) / h2
      kernel <- inv_sqrt_2pi / h2 * exp(-0.5 * diff^2)
      NijT2[i, j] <- NijT2[i, j] + kernel
    }
  }
  
  # 计算f1
  for (k in 1:p) {
    total <- 0
    for (i in 1:n) {
      for (j in 1:n) {
        if (i != j) {
          total <- total + Zij[i, j, k] * NijT2[i, j]
        }
      }
    }
    f1[k] <- total
  }
  
  # 牛顿迭代
  while (conv > 0.0001 && iter < max_iter) {
    iter <- iter + 1
    Eij <- matrix(0, n, n)
    
    # 计算Eij
    for (i in 1:n) {
      for (j in 1:n) {
        if (i != j) {
          temp <- sum(Zij[i, j, ] * xk)
          Eij[i, j] <- exp(temp)
        }
      }
    }
    
    # 计算f
    f <- numeric(p)
    for (k in 1:p) {
      total <- 0
      for (i in 1:n) {
        for (j in 1:n) {
          if (i != j) {
            total <- total + Zij[i, j, k] * Eij[i, j]
          }
        }
      }
      f[k] <- -total + f1[k]
    }
    
    # 计算D矩阵
    D <- matrix(0, p, p)
    for (k1 in 1:p) {
      for (k2 in 1:p) {
        total <- 0
        for (i in 1:n) {
          for (j in 1:n) {
            if (i != j) {
              total <- total + Zij[i, j, k1] * Zij[i, j, k2] * Eij[i, j]
            }
          }
        }
        D[k1, k2] <- total
      }
    }
    
    # 解线性方程组
    tryCatch({
      dd <- solve(D, f)
    }, error = function(e) {
      stop("Matrix D is singular or ill-conditioned.")
    })
    
    # 更新估计值
    xk <- xk + 0.05 * dd
    conv <- norm(as.matrix(dd), "F")  # Frobenius范数
    # cat(sprintf("Iteration %d: conv = %f\n", iter, conv))
    # print(round(xk, 4))
  }
  
  if (iter >= max_iter) {
    warning("Maximum iterations reached without convergence.")
  }
  
  return(xk)
}
