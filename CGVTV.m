%% The preprocess using vectorial total variation
function [ xMat,Iter ] = CGVTV( inMat,mask,tol,maxIter,lambda,noised )

[width,height,channel] = size(inMat);

Scale = max(max(max(inMat)));
inMat = inMat/Scale;
yMat = inMat.*mask;

% initialization
xMat = zeros(size(inMat));

grad = grad_2norm(mask,xMat,yMat)+lambda*grad_tvnorm(xMat);
dire = -grad; %search direction

alpha0 = 1;
beta = 0.6;
tau = 0.1;
k = 1;
epsilon = 1e-6;
nscale = 0.9;
    
% error1(k) = norm(grad,'fro');
error2 = (norm(grad(:,:,1),'fro')...
    + norm(grad(:,:,2),'fro')...
    + norm(grad(:,:,3),'fro'))...
    / norm(xMat(:,:,1)+epsilon,'fro')...
    + norm(xMat(:,:,2)+epsilon,'fro')...
    + norm(xMat(:,:,3)+epsilon,'fro');
    
% F_new = 0;
% while(error2(k)>tol && k<maxIter)
while(error2>tol && k<=maxIter)
    alpha = alpha0;
    num = 0;
    
    % the Armijo line search
    while(f_2norm(mask,xMat+alpha*dire,yMat)+...
            lambda*f_tvnorm(xMat+alpha*dire))>...
            (f_2norm(mask,xMat,yMat)+...
            lambda*f_tvnorm(xMat)+tau*alpha*(grad.*dire))
            %lamda*f_tvnorm(xMat)+tau*alpha*real(conj(grad).*dire))        
        alpha = beta*alpha;
        num = num+1;
    end
    if num>2
        alpha0 = beta*alpha0;
    end
    if num<1
        alpha0 = alpha0/beta;
    end
    
    xMat_old = xMat;
    xMat = xMat + alpha*dire;
    %xMat = max(xMat + alpha*dire,0);
    grad_old = grad;
    
    % the computation of the next search direction using conjugate gradient
    grad = grad_2norm(mask,xMat,yMat)+lambda*grad_tvnorm(xMat);
    new_norm = 0;
    old_norm = 0;
    for i = 1:channel
        new_norm = new_norm + norm(grad(:,:,i),'fro')^2;
        old_norm = old_norm + norm(grad_old(:,:,i),'fro')^2;
    end
    gamma = new_norm/old_norm;
%         gamma = gamma + norm(grad(:,:,i),'fro')^2/norm(grad_old(:,:,'fro')^2;
    dire = -grad+gamma*dire;
    
    k = k+1;
       
    error2 = (norm(grad(:,:,1),'fro')...
        + norm(grad(:,:,2),'fro')...
        + norm(grad(:,:,3),'fro'))...
        / (norm(xMat_old(:,:,1),'fro')...
        + norm(xMat_old(:,:,2),'fro')...
        + norm(xMat_old(:,:,3),'fro'));
%     error2
    lambda = lambda*nscale;
end
Iter = k - 1;
xMat = Scale*xMat;

end
    
        
%% the current value of term 1/2*||Ax-y||_2^2
function f_value = f_2norm(mask,X,Y)
f_value = 0;
for i = 1:size(X,3)
    f_value = f_value + norm(mask(:,:,i).*X(:,:,i)-Y(:,:,i),'fro')^2;
end
end

%% the current value of tv norm
function f_value = f_tvnorm(X)

X = [X X(:,end,:)];
X = [X;X(end,:,:)];
dx = (X(1:end-1,2:end,:)-X(1:end-1,1:end-1,:));
dy = (X(2:end,1:end-1,:)-X(1:end-1,1:end-1,:));

f_value = sum(sum(sum(sqrt(dx.^2+dy.^2))));
end

%% the gradient of the 1/2*||Ax-y||_2^2
function f_grad = grad_2norm(mask,X,Y)
f_grad = 2*mask.*(mask.*X-Y);
end

%% the gradient of tv norm
function f_grad = grad_tvnorm(X)

s = size(X);
delta = 1e-14;

X = [X(:,1,:) X X(:,end,:)];
X = [X(1,:,:);X;X(end,:,:)];

% the vertical and horizontal gradient of point (i,j)
h1 = [X(2:end-1,3:end,:)-X(2:end-1,2:end-1,:)];
v1 = [X(3:end,2:end-1,:)-X(2:end-1,2:end-1,:)];

% the vertical and horizontal gradient of point(i-1,j) (the above point)
h2 = [X(1:end-2,3:end,:)-X(1:end-2,2:end-1,:)];
v2 = [X(2:end-1,2:end-1,:)-X(1:end-2,2:end-1,:)];

% the vertical and horizontal gradient of point(i,j-1) (the left point)
h3 = [X(2:end-1,2:end-1,:)-X(2:end-1,1:end-2,:)];
v3 = [X(3:end,1:end-2,:)-X(2:end-1,1:end-2,:)];

grad1 = sqrt(sum(h1.^2 + v1.^2,3) + delta);
grad2 = sqrt(sum(h2.^2 + v2.^2,3) + delta);
grad3 = sqrt(sum(h3.^2 + v3.^2,3) + delta);

f_grad = zeros(s);
for c = 1:size(X,3)
    f_grad(:,:,c) = -v1(:,:,c)./grad1 - h1(:,:,c)./grad1 + v2(:,:,c)./grad2 + h3(:,:,c)./grad3;
end
end