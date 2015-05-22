#' @title Unrestricted, unstandardized covariance-based PCA with Varimax 
#' rotation for Event-Related Potentials (ERPs).
#' 
#' @description This package implements the set of recommendations laid out by
#' Jurgen Kayser and Craig E. Tenke on a series of publications (see references 
#' below) on how to optimize the decomposition of Event-Related Potential (ERP)
#' data into useful latent factors. These authors propose a PCA decomposition of
#' an ERP data matrix based on its unstandardized covariance matrix, followed by
#' a Varimax rotation of all of its principal components, in an attempt to find 
#' simpler latent factors that may correspond to the event-related potentials 
#' that researchers using EEG/MEG are interesting in measuring.
#' 
#' The proposed workflow involves the decomposition of the ERP matrix in its 
#' principal components using the unstandardized covariance matrices of the 
#' original ERP data matrix followed by a Varimax rotation of all of the 
#' principal components (what Kayser and Tenke call an ``unrestricted'' 
#' solution).
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
#' @source This package is a port of original MATLAB functions created by Jurgen
#'   Kayser, which can be freely downloaded from 
#'   \url{http://psychophysiology.cpmc.columbia.edu/mmedia/Kayser2003a/Appendix.html}
#'   
#'
#' @name erpPCA-package
#' @docType package
#' @author Diogo Almeida \email{diogo@@nyu.edu}
#' @keywords package
NULL
