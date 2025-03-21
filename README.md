This repository contains the source code for the golem95 library. 
It can perform the tensor reduction of one-loop integrals for up to six external legs and also contains the one-loop master integrals. 
For integrals with massless propagators, it contains a method to evaluate tensor integrals directly via a fast one-dimensional numerical integration, this way avoiding the occurrence of inverse Gram determinants. This improves the numerical stability a lot in kinematic regions where the Gram determinant tends to zero.
