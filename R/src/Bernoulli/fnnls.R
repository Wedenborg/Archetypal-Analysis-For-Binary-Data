fastnnls <- function(X, Xty, tol, b, PP, device = "cpu") {
  # Number of columns in the input data X.
  n <- ncol(X)
  
  if (length(PP) == 0) {
    P <- rep(NA, n)
    PP <- which(!is.na(P))
    
    Z <- 1:n # Initialize the set of candidate indices.
    ZZ <- which(is.na(P))
    t = rep(0,nrow(X))
  } else {
    P <- rep(NA, n)
    P[PP] <- PP
    
    Z <- 1:n # Initialize the set of candidate indices.
    Z[PP] <- NA
    ZZ <- which(is.na(P))
    if (length(PP)==1){
      t <- X[, PP] * x[PP]
      
    } else{
      t <- X[, PP] %*% x[PP]
      
    }
  }
  
  x <- as.double(b)
  z <- as.double(b)
  
  ij <- NULL
  iterOuter <- 0
  
  scale <- mean(colSums(X * X))
  lambda1 <- 1e9 * scale
  lambda2 <- 1e-6 * scale
  #t <- X[, PP] %*% x[PP]
  w <- Xty + lambda1 - t(X) %*% t - sum(x) * lambda1 - lambda2 * x
  iter <- 0
  
  if (length(PP) == 0) {
    XtX_PP <- NULL
  } else {
    XtX_PP <- t(X[, PP]) %*% X[, PP]
  }
  
  while (length(ZZ) > 0 && sum(w[ZZ] > tol) > 0 && iterOuter < 1000) {
    iter <- 0
    temp <- w[ZZ]
    t <- which.max(temp)
    t <- ZZ[t]
    
    PP <- c(PP, t)
    ZZ <- ZZ[ZZ != t]
    
    if (is.null(XtX_PP)) {
      XtX_PP <- t(X[, PP]) %*% X[, PP]
    } else {
      if (length(PP[-length(PP)]) == 1) {
        XtX_PP <- cbind(XtX_PP[1], t(X[, t]) %*% X[, PP[-length(PP)]])
        XtX_PP <- rbind(XtX_PP, t(t(X[, PP]) %*% X[, t]))
      } else {
        XtX_PP <- cbind(XtX_PP, t(t(X[, t]) %*% X[, PP[-length(PP)]]))
        XtX_PP <- rbind(XtX_PP, t(t(X[, PP]) %*% X[, t]))
      }
    }
    z <- rep(0, length(z))
    z[PP] <- solve(XtX_PP + lambda1 + diag(length(PP)) * lambda2, Xty[PP] + lambda1)
    
    while (sum(z[PP] <= tol) > 0 && iter < 500) {
      temp <- which(z[PP] <= tol)
      QQ <- PP[temp]
      
      a <- x[QQ] / (x[QQ] - z[QQ])
      a[is.na(a)] <- Inf
      alpha <- min(a)
      
      x <- x + alpha * (z - x)
      
      t1 <- which(x[PP] < tol)
      ij <- PP[t1]
      
      idx_PP <- if (length(ij) == 0) rep(TRUE, length(PP)) else !PP %in% ij
      PP <- PP[idx_PP]
      
      ZZ <- c(ZZ, ij)
      
      if (is.null(XtX_PP)) {
        XtX_PP <- t(X[, PP]) %*% X[, PP]
      } else {
        XtX_PP <- XtX_PP[idx_PP, idx_PP]
      }
      
      z <- rep(0, length(z))
      z[PP] <- solve(XtX_PP + lambda1 + diag(length(PP)) * lambda2, Xty[PP] + lambda1)
      
      iter <- iter + 1
    }
    
    x <- z
    t <- if (length(PP) == 1)  X[, PP] * x[PP] else  X[, PP] %*% x[PP]
    w <- Xty + lambda1 - t(X) %*% t - sum(x) * lambda1 - lambda2 * x
    
    iterOuter <- iterOuter + 1
  }
  
  out <- z
  #return(list(out = out, w = w, PP = PP))
  return(out)
}
