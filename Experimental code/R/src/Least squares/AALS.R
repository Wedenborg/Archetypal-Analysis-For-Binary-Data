
AALS <- function(X, n_samples, n_arc, C = NULL, S = NULL, gridS = FALSE) {
  if (is.null(C)) {
    # Generate random numbers
    C <- abs(log(matrix(runif(n_arc * n_samples), nrow = n_samples, ncol = n_arc)))
    
    # Normalize the matrix
    C <- t(t(C) / colSums(C))
  }
  
  if (is.null(S)) {
    # Generate random numbers
    S <- abs(log(matrix(runif(n_arc * n_samples), nrow = n_arc, ncol = n_samples)))
    # Normalize the matrix
    S <- t(t(S) / colSums(S))
  }
  
  L <- numeric()
  EV <- numeric()
  iter <- 0
  loss <- 0
  loss_old <- 1e10
  SSt <- S %*% t(S)
  
  SST <- sum(X^2)
  
  XC <- X %*% C
  
  CtXtXC <- t(XC) %*% XC
  XCtX <- t(XC) %*% X
  
  PPP <- vector("list", n_arc)
  LS = 0
  LC = 0 
  
  while (abs(loss_old - loss) > 1e-6 * abs(loss_old) && iter < 10) {
    loss_old <- loss
    
    S_update <- S_updateR(S, CtXtXC, XCtX, SSt, n_arc, n_samples, gridS)
    S <- S_update$S
    SSt <- S_update$SSt
    Ls = S_update$L
    
    LS = c(LS,Ls)
    
    C_update <- C_updateR(X, C, S, SSt, XC, CtXtXC, n_arc, PPP)
    C <- C_update$C
    XC <- C_update$XC
    CtXtXC <- C_update$CtXtXC
    PPP <- C_update$PP
    
    Lc = C_update$L
    LC = c(LC,Lc)
    
    
    XCtX <- t(XC) %*% X
    
    loss <- SST - 2 * sum(XCtX * S) + sum(CtXtXC * SSt)
    
    L <- c(L, loss)
    
    varExpl <- (SST - loss) / SST
    
    EV <- c(EV, varExpl)
    
    iter <- iter + 1
  }
  
  return(list(C = C, S = S, L = L, EV = EV, LS = LS, LC = LC ))
}
