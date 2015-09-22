% This is a re-implemented Matlab source code of paper:
% Wei Li, Lei Zhao, Duanqing Xu, and Dongming Lu, "A Non-local Method for Robust Noisy Image 
% Completion", ECCV 2014(IV:61-74).
% The author of the code is Wei Li, 05/07/2015.
% Contact: leewei.david@gmail.com

clear;clc; close all;
oImg = double(imread('lena.png'));%original image
[height width channel] = size(oImg);

sample_rate = 0.70;
sample_num = round(height*width*sample_rate);
mask = randsample(height*width,sample_num);
mask_matrix = zeros(height,width);
mask_matrix(mask) = 1;

noise_std = 25;
nImg = oImg +  noise_std*randn(height,width,channel);%noised image

I = nImg;
mask3 = zeros(height,width,channel);
for i = 1:channel
    I(:,:,i) = nImg(:,:,i).*mask_matrix;
    mask3(:,:,i) = mask_matrix;
end

%% preprocessing using robust total variation completion method with noised input
maxIter = 500;
tol = 1e-4;
lambda = 0.2;
disp('preprocessing...');
% Treat the input as noisy image or not
[pre_result,n] = CGVTV(I,mask3,tol,maxIter,lambda,0); 

%% searching for k nearest neighbour patches
cores = 2;
source = uint8(pre_result);
patch_size = 7;
patch_num = 20;
nn_iters = 50;
disp('executing patchmatch...');
% the original version of PatchMatch
knn = nnmex(source, source, 'gpucpu', patch_size, nn_iters, [], [], [], [], cores, [], [], [], [], [], patch_num);

%% forming low-rank matrix completion problem
sum_result = zeros(height,width,channel);
count = zeros(height,width,channel);
c = 0;
% the boundry process is ignored for simplicity
for y = 1:3:height-patch_size+1
    for x = 1:3:width-patch_size+1
        
        lrm_h = patch_size*patch_size*channel; % the height of the to be formed low rank matrix
        lrm_w = patch_num;
        lrm = zeros(lrm_h,lrm_w);
        lrm_mask = zeros(lrm_h,lrm_w);
        current = find(knn(y,x,3,:) == min(knn(y,x,3,:)));
        
        % cancatenate each patch into a vector and stack all the patches
        for n = 1:lrm_w
            tx = knn(y,x,1,n)+1;
            ty = knn(y,x,2,n)+1;
            lrm(:,n) = reshape(I(ty:ty+patch_size-1,tx:tx+patch_size-1,:),lrm_h,1);
            lrm_mask(:,n) = reshape(mask3(ty:ty+patch_size-1,tx:tx+patch_size-1,:),lrm_h,1);
        end

        mask = lrm_mask;
        sigma = noise_std;
        initial = [];
        
        % performance largely depends on these parameters, which are not optimized and need
        % to be studied in the future
        mr = sum(sum(mask == 0))/(lrm_h*lrm_w); %missing rate
        beta = 0.45/(sqrt(lrm_h*lrm_w)*mr); 
        gamma = 1.6;
        tol2 = 1e-5;
        delta = sqrt(sum(sum(mask))) * sigma * 0.8;
        [lrm_result,iter,error] = ADM_noise(lrm,mask,initial,beta,gamma,tol2,delta);
     
        sum_result(y:y+patch_size-1,x:x+patch_size-1,:) = sum_result(y:y+patch_size-1,x:x+patch_size-1,:)...
            + reshape(lrm_result(:,current(1)),patch_size,patch_size,channel);
        count(y:y+patch_size-1,x:x+patch_size-1,:) = count(y:y+patch_size-1,x:x+patch_size-1,:)...
            + ones(patch_size,patch_size,channel);

    end
end

% the result here tends to be a little darker, so a small bias of the color is
% needed for the purpose of display in the paper. 
% The reason is still ambiguous now.
rImg = sum_result./count + 5;
figure;
imshow(rImg/255);
title('final result');

