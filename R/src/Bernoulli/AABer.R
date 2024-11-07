

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
  lambda2 <- 1e-4 * scale
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

# LossFunc function
LossFunc <- function(X, PC, S) {
  loss <- sum(-X * log(PC %*% S)) - sum((1 - X) * log(1 - (PC %*% S)))

  return(loss)
}


SUpdate <- function(X, P, PC, S, SSt, n_arc, n_samples) {
  alpha <- rep(0, n_samples)
  tol <- 1e-12

  # Create a grid of all possible combinations of archetypes
  # Create a grid of all possible pairs of archetypes
  grid <- combn(n_arc, 2)
  grid <- grid[, sample(ncol(grid))]

  lS = LossFunc(X,PC,S)


  for (j in seq(ncol(grid))) {

    SS <- colSums(S[grid[, j], ])
    S0 <- S[grid[, j], ]

    K <- -t(PC[, grid[, j]]) %*% (X / (PC %*% S)) + t(PC[, grid[, j]]) %*% ((1 - X) / (1 - (PC %*% S)))

    QuadS <- einsum::einsum('ji,jk,jl->ikl',PC[, grid[, j]],PC[, grid[, j]],X / ((PC %*% S) ^ 2) + (1 - X) / ((1 - (PC %*% S)) ^ 2))

    linS <- K - 2*einsum::einsum('kl,ikl->il',S0,QuadS)


    denominator <- 2 * SS[SS > tol] * (QuadS[1, 1, ] - QuadS[1, 2, ] - QuadS[2, 1, ] + QuadS[2, 2, ])[SS > tol]
    nominator <- SS[SS > tol] * (QuadS[1, 2, ] + QuadS[2, 1, ] - 2 * QuadS[2, 2, ])[SS > tol] + linS[1, SS > tol] - linS[2, SS > tol]

    # Compute the alpha values.
    alpha[SS > tol] <- - (nominator / denominator)
    alpha[alpha < 0] <- 0
    alpha[alpha > 1] <- 1

    S[grid[1, j], SS > tol] <- SS[SS > tol] * alpha[SS > tol]
    S[grid[2, j], SS > tol] <- SS[SS > tol] * (1 - alpha[SS > tol])

    SSt[grid[, j], ] <- S[grid[, j], ] %*% t(S)


    SSt[, grid[, j]] <- t(SSt[grid[, j], ])
    lS = c(lS, LossFunc(X,PC,S))
  }

  return(list(S = S, SSt = SSt, lS = lS))
}



CUpdate <- function(X, P, C, S, PC, n_arc, PPP) {
  for (i in seq(n_arc)) {
    v <- as.vector(((X / ((PC %*% S) ^ 2)+1e-6) %*% (S[i, ] ^ 2)) + (((1 - X) / ((1 - (PC %*% S)) ^ 2)+1e-6) %*% (S[i, ] ^ 2)))
    grad <- as.matrix(t(P) %*% ((X / (PC %*% S)+1e-6) %*% S[i, ]) - t(P) %*% (((1 - X) / (1 - PC %*% S)+1e-6) %*% S[i, ]))

    xtilde <- P*sqrt(v)


    Xty <- grad + t(xtilde) %*% (xtilde %*% C[, i])
    tol <- 1e-9
    C_k <- fastnnls(xtilde, Xty, tol, C[, i], c())

    C[, i] <- C_k
    PC[, i] <- P %*% C_k
  }

  lC = LossFunc(X,PC,S)
  return(list(C = C, PC = PC, PPP = PPP, lC = lC))
}


# ClosedFormArchetypalAnalysis function
#' @export
ClosedFormArchetypalAnalysis <- function(X, n_arc, maxIter = 500, convCriteria = 1e-6, C = NULL, S = NULL, init = FALSE) {
  eps <- 1e-3
  n_features <- nrow(X)
  n_samples <- ncol(X)

  P <- X + eps - 2 * eps * X

  iter <- 0
  loss <- 9
  lossOld <- 10
  L <- c()

  if (is.null(C)) {
    # Generate random numbers
    C <- log(matrix(runif(n_arc * n_samples), nrow = n_samples, ncol = n_arc))

    # Normalize the matrix
    C <- t(t(C) / colSums(C))

  }

  if (is.null(S)) {
    # Generate random numbers
    S <- log(matrix(runif(n_arc * n_samples), nrow = n_arc, ncol = n_samples))
    # Normalize the matrix
    S <- t(t(S) / colSums(S))

  }

  if (init) {
    C <- AALS(X, n_samples, n_arc, C, S, gridS = TRUE, maxIter = 50)[[1]]
    S <- AALS(X, n_samples, n_arc, C, S, gridS = TRUE, maxIter = 50)[[2]]
  }

  PC <- P %*% C
  SSt <- S %*% t(S)
  PPP <- vector("list", n_arc)

  while (abs(loss - lossOld) / abs(loss) > convCriteria && iter < 1000) {
    lossOld <- loss

    # Update S
    S_update <- SUpdate(X, P, PC, S, SSt, n_arc, n_samples)
    S <- S_update$S
    SSt <- S_update$SSt

    # Update C
    C_update <- CUpdate(X, P, C, S, PC, n_arc, PPP)
    C <- C_update$C
    PC <- C_update$PC
    PPP <- C_update$PPP

    # Update loss
    loss <- LossFunc(X, PC, S)
    L <- c(L, loss)
    iter <- iter + 1
  }

  return(list(C = C, S = S, L = L))
}
