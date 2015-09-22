function [ output,iter,error ] = ADM_noise( input,mask,initial,beta,gamma,tol, delta )
%ADM_NOISE Summary of this function goes here
%   Detailed explanation goes here

[m,n] = size(input);
%% initilization
maxIter = 500;
mask = (mask == 0);%here 1 stands for empty
X = input;
Y = X;
if isempty(initial)
    X(mask) = 0;
    Y(mask) = 0;
else
    X(mask) = initial(mask);
    Y(mask) = initial(mask);
end
Z = X;
%% iteration
for i = 1:maxIter
    %calculate y
    B = X - Z/beta;
    T = B - input;
    T(mask) = 0;
   % norm(T,'fro')
   % a = min(sigma/norm(T,'fro'),1);
    Y = (min(delta/norm(T,'fro'),1)-1)*T+B;
    %Y(mask) = B(mask);
    
    %calculate x
    A = Y + Z/beta;
    newX = findX(A,beta);
    
    error(i) = norm(newX-X,'fro')/norm(X,'fro');
    if error(i) < tol
        break;
    end
    X = newX;
    
    %calculate z
    Z = Z - gamma*beta*(X-Y);
end
output = X;
iter = i;
end


