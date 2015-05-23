# Version 1.0 (May 2015)

* Full port of erpPCA.m and related MATLAB functions from Jurgen Kayser completed.

* There are a couple of options in the package that go beyond the simple porting of the original MATLAB functions. These include
    * Decompose the data matrix with eigen or svd (original functions use eigen)
    * Two options for estimating the rank of the svd matrix (original, which uses svd of correlation matrix, and qr, which does not call svd, and therefore is probably more computationally efficient.). The options are not interchangeable, though.
    * Two options for Varimax rotation. The original solution uses the dedicated function DoVarimax4M, but the user may also use the base R varimax function.