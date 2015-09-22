function [L_hat,S_hat] = rPCA(patchdata,psite,sigma)
% patchdata is the noised incomplete low rank matrix
% psite is the position set of all the empty elements
% sigma is the standard deviation of noise


[nr,nc] = size(patchdata);
p = nnz(psite)/(nr*nc);
%lambda = 1*1/sqrt(max(nr,nc));
lambda = 0.8*1/sqrt(max(nr,nc));
mu = 1*sqrt(p)*sigma*(sqrt(nr)+sqrt(nc));
% mu = 2*(sqrt(nc)+sqrt(nc))*sigma;
maxIter = 30;
tol = 1e-2;
[L_hat,S_hat,iter] = partial_proximal_gradient_rpca2(patchdata, psite, lambda, mu, maxIter, tol);

end