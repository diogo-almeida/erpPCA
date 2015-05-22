#' @title Performs unrestricted, unstandardized covariance-based PCA followed by
#'   Varimax rotation
#' 
#' @description This function implements a PCA decomposition of a data matrix 
#' based on its unstandardized covariance matrix, followed by a Varimax rotation
#' of all of its principal components, in an attempt to find simpler latent 
#' factors.
#' 
#' @details This workflow of matrix decomposition followed by Varimax rotation 
#'   has been advocated for electrophysiological studies on humans using the 
#'   Event-Related Potential (ERP) technique by Jürgen Kayser and Craig E. 
#'   Tenke. These authors suggest that very good results can be obtained by
#'   using (1) unstandardized covariance matrices of the original data matrix
#'   followed by (2) a Varimax rotation of all of the principal components (what
#'   they call an ``unrestricted'' rotation). The arguments are laid out in a
#'   series of papers (see references below)
#'  
#' This function is a R port of the original MATLAB function erpPCA.m created
#' by Jürgen Kayser. According to the original documentation, it implements
#' ``the PCA agorithms used by BMDP-4M (Dixon, 1992) and SPSS 10.0 FACTOR''. The
#' original MATLAB code can be freely downloaded at 
#' \url{http://psychophysiology.cpmc.columbia.edu/mmedia/Kayser2003a/Appendix.html}
#' 
#' @param X A matrix. This is the data matrix to be decomposed into its
#'   principal components. For a temporal decomposition (ie., latent factors
#'   will correspond to time courses), the matrix should be organized
#'   with ERP waveforms of each sensor in rows (cases) and time samples in 
#'   columns (variables). For spatial decomposition (ie., latent factors
#'   will correspond to scalp topographies), the matrix should be organized
#'   with time samples in rows (cases) and sensors in columns (variables).
#'   
#' @return This function returns a list containing the unrotated principal
#'   component loadings, the Varimax rotated factor loadings in a 
#'   variables-by-factors matrix, the factor scores obtained for the rotated
#'   factors in a cases-by-factors matrix. It also outputs the Eingenvalues
#'   of the covariance matrix and the their explained variance before and after
#'   rotation.
#' 
#' \item{Unrotated}{Matrix containing the unrotated principal components
#' loadings.}
#' 
#' \item{Rotated}{Matrix containing the varimax rotated factor loadings}
#' 
#' \item{Factor.Scores}{Matrix containing the factor scores of the varimax
#' rotated factor loadings.}
#' 
#' \item{Variance}{Matrix containing the Eigenvalues and explained variance
#' before and afte rotation.}
#'   
#' @references Kayser, J., & Tenke, C. E. (2003). Optimizing PCA methodology for
#'   ERP component identification and measurement: theoretical rationale and 
#'   empirical evaluation. Clinical neurophysiology, 114(12), 2307-2325.
#'   
#' @references Kayser, J., & Tenke, C. E. (2005). Trusting in or breaking with 
#'   convention: towards a renaissance of principal components analysis in 
#'   electrophysiology. Clinical Neurophysiology, 116(8), 1747-1753.
#'   
#' @references Kayser, J., & Tenke, C. E. (2006). Consensus on PCA for ERP data,
#'   and sensibility of unrestricted solutions. Clinical Neurophysiology, 
#'   117(3), 703-707.
#'   
#' @seealso \code{\link[stats]{princomp}}, \code{\link[stats]{prcomp}}, 
#'   \code{\link[stats]{factanal}} from package \code{stats} for other basic 
#'   implementations of principal component analysis and factor analysis.
#'   
#' @source
#'   \url{http://psychophysiology.cpmc.columbia.edu/mmedia/Kayser2003a/Appendix.html}
#'   
#' @export
#'    
#' @examples
#' ## using the iris dataset and principal components analysis
#' erpPCA(iris[, 2:4])
#'
erpPCA <- function(X) {
##############
# Port the erpPCA.m file from Kayser & Tenke 2003 to R
# Diogo Almeida - 6/21/12
##############

  # Helper function definition
  # ------------------------------
  # repmat (emulate repmat from matlab)
  # ------------------------------
  repmat <- function(a, n, m) {
    kronecker(matrix(1, n, m), a)
  }

  # Input checks
  if (is.list(X)) {
    cat("Your data needs to be input as a matrix. You input a list.\n")
    cat("I'm converting your data to the appropriate format now.\n")
    X <- as.matrix(X)
    if (is.double(X)) {
      cat("Your data is now a matrix!\n")
    }
  }
    
  D <- cov(X)
  eigenD <- eigen(D)
  EM <- eigenD$vectors
  EV <- eigenD$values
  UL <- EM %*% sqrt(diag(EV))                  # svd(X)$v gives different values
  rk <- sum(svd(cor(X))$d > 1e-4)
  u <- sort(EV, decreasing = T)[1:rk]
  LU <- UL[, 1:rk]
  s <- matrix(1, ncol = rk)
  s[abs(apply(LU, 2, max)) < abs(apply(LU, 2, min))] <- -1
  LU <- LU * repmat(s, dim(LU)[1], 1)
  RL <- DoVarimax4M(LU)$Y
  EVr <- colSums(RL * RL)
  r   <- sort(EVr, decreasing = T)
  rx  <- order(EVr, decreasing = T)
  LR  <-  RL[, rx]
  s <- matrix(1, ncol =  ncol(LR))
  s[abs(apply(LR, 2, max)) < abs(apply(LR, 2, min))] <- -1
  LR <- LR * repmat(s, ncol(X), 1) 
  tv <- sum(EV)
  VT <- matrix(c(u, (100 * u) / tv, r, (100 * r) / tv), byrow = F, 
               nrow = length(u))
  FSCFr <- LR %*% solve(t(LR) %*% LR)
  FSCFr <- FSCFr * repmat(sqrt(diag(D)), 1, rk) 
  mu <- apply(X, 2, mean)
  sigma <- apply(X, 2, sd)
  Xc <- sweep(X, 2, mu)
  FSr = matrix(0, nrow = nrow(X), ncol = ncol(X))
  for (n in 1:nrow(X)) {
    for (m in 1:rk) {
      FSr[n, m] <- sum((Xc[n, ] / sigma) * FSCFr[, m])
    }
  }
  
  # Prettify output (VT)
  factor.labels <- paste("Factor", 1:nrow(VT), sep = ".")
  header.labels <- c("Eigenvalue.Before", "Var.Explained.Before", 
                     "Eigenvalue.After", "Var.Explained.After")
  rownames(VT) <- factor.labels
  colnames(VT) <- header.labels
  
  # return the same output
  output <- list(Unrotated = LU, Rotated = LR, Factor.Scores = FSr, 
                 Variance = VT)
  output
}
