clc;
clear all;
I = imread('image.jpg');
%imshow(I);
I_G = rgb2gray(I);
%imshow(I_G);
[M N] = size(I_G);
a = mod(M,8);
b = mod(N,8);
MM = M - a;
NN = N - b;
I_p = double(I_G);
I_p = I_G(1:MM,1:NN);
I_re = zeros(MM,NN);

    
for i = 1:8:MM
    for j = 1:8:NN
        K = I_p(i:i+7,j:j+7);
        K_DCT = dct2(K);
        K_DCT(abs(K_DCT)<(100)) = 0;
        
              
        K_RE = idct2(K_DCT);
        I_re(i:i+7,j:j+7) = K_RE;
        
    end
end
 imshow(uint8(I_re));
 figure,imshow(I_p);
 figure, imshow(uint8(I_re)- I_p);
% err = (uint8(I_re) - I_p);
% err = err.*err;
% MSE(k+1) = sum(sum(err))/(MM*NN)
% Comp_ratio(k+1) = ((MM*NN) - no_of_zeroed)/(MM*NN);


% plot(MSE,Comp_ratio);
