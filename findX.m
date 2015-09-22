function [ X,rank ] = findX( A,beta )
%singular value shrinkage
%   Detailed explanation goes here
% [S,V,D] = svd(A);
% n = sum(diag(V) > 1/beta);
% %n = sum(diag(V) > beta);
% S = S(:,1:n);
% V = V(1:n,1:n);
% D = D(:,1:n);
% 
% %V = max(diag(V)-beta,0);
% V = max(diag(V)-1/beta,0);
% r = sum(V>0);
% X = S(:,1:r) * diag(V(1:r)) * D(:,1:r)';
% end
%% test for another processing of find x
[S,V,D] = svd(A);
rank = sum(diag(V)>0);

%V = max(diag(V)-beta,0);
V = max(diag(V)-1/beta,0);
r = sum(V>0);
X = S(:,1:r) * diag(V(1:r)) * D(:,1:r)';