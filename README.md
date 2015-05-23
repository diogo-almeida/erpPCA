# erpPCA

R port of MATLAB functions erpPCA.m and others created by Jürgen Kayser. The original MATLAB functions can be freely downloaded from http://psychophysiology.cpmc.columbia.edu/mmedia/Kayser2003a/Appendix.html

## Scope

This package implements only the PCA + Varimax rotation workflow advocated by Kayser and Tenke (2003, 2005, 2006) for Event-Related Potential (ERP) data, and is not a general approach to PCA nor to factor analysis.

A more flexible factor analysis approach to ERP data can be found in the ERP PCA Toolkit by Joseph Dien (Dien, 2010a), which is available only in MATLAB. Perhaps a future version of this package will implement some aspects of this more flexible approach, namely:

* The ability to use different rotations (especially Promax, as recommended by Dien (2010b, Dien et al., 2007)

* The ability to extract and rotate only a user-defined number of factors.

The R ecosystem already provides plenty of PCA and factor analyses packages that can probably be co-opted for these purposes.

## References

* Dien, J., Khoe, W., & Mangun, G. R. (2007). Evaluation of PCA and ICA of simulated ERPs: Promax vs. Infomax rotations. Human brain mapping, 28(8), 742-763.

* Dien, J. (2010a). The ERP PCA Toolkit: An open source program for advanced statistical analysis of event-related potential data. Journal of neuroscience methods, 187(1), 138-145.

* Dien, J. (2010b). Evaluating two‐step PCA of ERP data with geomin, infomax, oblimin, promax, and varimax rotations. Psychophysiology, 47(1), 170-183.

* Kayser, J., & Tenke, C. E. (2003). Optimizing PCA methodology for ERP component identification and measurement: theoretical rationale and empirical evaluation. Clinical Neurophysiology, 114(12), 2307-2325.

* Kayser, J., & Tenke, C. E. (2005). Trusting in or breaking with convention: towards a renaissance of principal components analysis in electrophysiology. Clinical Neurophysiology, 116(8), 1747-1753.

* Kayser, J., & Tenke, C. E. (2006). Consensus on PCA for ERP data, and sensibility of unrestricted solutions. Clinical Neurophysiology, 117(3), 703-707.

