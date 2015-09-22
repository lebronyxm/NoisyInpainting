This is a re-implemented and simplified Matlab code of the paper:
Wei Li, Lei Zhao, Duanqing Xu, and Dongming Lu, "A Non-local Method for Robust Noisy Image Completion", ECCV 2014 (IV:61-74).
The author of the code is Wei Li, 05/08/2015.

The file list:
Demo.m        The entrance of the algorithm.
CGVTV.m       The pre-process of the algorithm using vectorial total variation (VTV).
ADM_noise.m   The solver (Alternating Direction Method of Multipliers) for the 
              optimization problem in the algorithm.
findX.m       An auxiliary function to ADMM
rPCA.m        A different solver for the optimization problem
nnmex.mexw64  The mex file of PatchMatch:
              http://gfx.cs.princeton.edu/gfx/pubs/Barnes_2009_PAR/index.php
              In this re-implemented version, the original PatchMatch is directly 
              borrowed rather than the improved version reported in the paper. 

This version is not exactly the same as the original one, and the parameters in the algorithm are not well tuned (a tricky process, which is considered to be the limitation of the algorithm). The main purpose of this implementation is to help readers to understand the processing of the algorithm more clearly. 

Contact of the author: leewei.david@gmail.com
Any reports about problems and bugs are welcome and grateful.