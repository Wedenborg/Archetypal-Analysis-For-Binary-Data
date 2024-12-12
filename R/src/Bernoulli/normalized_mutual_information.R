# Calculate mutual information
calcMI <- function(z1, z2) {
  eps <- 10e-16
  P <- z1 %*% t(z2)
  PXY <- P / sum(P)
  PXPY <- outer(rowSums(PXY), colSums(PXY))
  #ind <- which(PXY > 0)
  #MI <- sum(PXY[ind] * log(eps + PXY[ind] / (eps + PXPY[ind])))
  MI <- sum(PXY * log(eps + PXY / (eps + PXPY)))
  return(MI)
}

# Calculate normalized mutual information
calcNMI <- function(z1, z2) {
  NMI <- (2 * calcMI(z1, z2)) / (calcMI(z1, z1) + calcMI(z2, z2))
  return(NMI)
}
