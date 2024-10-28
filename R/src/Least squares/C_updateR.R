# C update function
C_updateR <- function(X, C, S, SSt, XC, CtXtXC, n_arc, PP) {
  loss <- list()
  L = 0
  
  for (i in 1:n_arc) {
    
    
    Xtilde <- X * sqrt(SSt[i, i])
    Xty <- 2 * (t(X) %*% (X %*% S[i, ]) - t(X) %*% (X %*% (C %*% (S %*% S[i, ])))) + t(Xtilde) %*% (Xtilde %*% C[, i])
    

    C_k_PP <- fastnnls(Xtilde,Xty,1e-9,C[, i],c())
    
    #C_k <- C_k_PP$x
    C_k <- C_k_PP
    #PP[i] <- 0 #PP[i]
    
    C[, i] <- C_k
    XC[, i] <- X %*% C_k
    
    CtXtXC[i, ] <- XC[, i] %*% XC
    CtXtXC[, i] <- CtXtXC[i, ]
    
    SST = sum(X^2)
    
    loss <- SST - 2 * sum(t(XC)%*%X * S) + sum(CtXtXC * SSt)
    
    L <- c(L, loss)
    
  }
  
  return(list(C = C, XC = XC, CtXtXC = CtXtXC, PP = PP,L = L))
}
