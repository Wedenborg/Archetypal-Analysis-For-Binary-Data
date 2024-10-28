# S update function
S_updateR <- function(S, CtXtXC, XCtX, SSt, n_arc, n_samples, gridS = FALSE) {
  tol <- 1e-3
  eta <- 1e-16
  
  # Create a grid of all possible pairs of archetypes
  grid <- combn(n_arc, 2)
  grid <- grid[, sample(ncol(grid))]
  
  
  if (gridS) {
    grid <- cbind(grid, grid)
  }
  
  alpha <- rep(0, n_samples)
  L = 0 
  for (j in 1:ncol(grid)) {
    SS <- colSums(S[grid[, j], ])
    S0 <- S[grid[, j], ]
    
    hess <- CtXtXC[grid[, j], grid[, j]]
    h_lin <- CtXtXC[grid[, j], ]
    
    
    d1 <- -2 * XCtX[grid[, j], ] + 2 * (h_lin %*% S)
    d2 <- 2 * (hess %*% S0)
    d <- d1 - d2
    denominator <- 2 * SS[SS > tol] * (hess[1, 1] - hess[1, 2] - hess[2, 1] + hess[2, 2])
    nominator <- SS[SS > tol] * (hess[1, 2] + hess[2, 1] - 2 * hess[2, 2]) + d[1, SS > tol] - d[2, SS > tol]
    
    # Compute the alpha values
    alpha[SS > tol] <- - t(t(nominator) / (denominator + eta))
    
    # Clip the values of alpha to the interval [0, 1]
    alpha[alpha < 0] <- 0
    alpha[alpha > 1] <- 1
    
    # Update the archetypes using the computed alpha values
    S[grid[1, j], SS > tol] <- SS[SS > tol] * alpha[SS > tol]
    S[grid[2, j], SS > tol] <- SS[SS > tol] * (1 - alpha[SS > tol])
    
    
    SSt[grid[, j], ] <- S[grid[, j], ] %*% t(S)
    SSt[, grid[, j]] <- t(SSt[grid[, j], ])
    
    SST = sum(X^2)
    
    loss <- SST - 2 * sum(XCtX * S) + sum(CtXtXC * SSt)
    
    L <- c(L, loss)
  }
  
  return(list(S = S, SSt = SSt,L = L))
}   
