clc;
clear all;
q=1;

%%%%%%%%%%%%%%%%%% Read frame %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for frmno=0:10;
  myfilename = sprintf('mother/%d.pgm',frmno);
  framR(:,:,q)=(imread(myfilename));
  q=q+1;
end

%%%%%%%%%%%%%%%%% Individuall frame dct %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:6;
    
    I_G=framR(:,:,k);
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
            I_re(i:i+7,j:j+7) = K_DCT;
         end
    end
    frmdct(:,:,k)=I_re;    
    frmdct(:,:,k)=dct2(I_p(:,:));
end

%%%%%%%%%%%%%%%%%%%%%%%  Delta calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for k=1:6
    
    I_G=framR(:,:,k);
    [M N] = size(I_G);
    a = mod(M,8);
    b = mod(N,8);
    MM = M - a;
    NN = N - b;
    I_p = double(I_G);
    I_p = I_G(1:MM,1:NN);
    delt=zeros(MM,NN);
    
    for i=1:MM
        for j=1:NN
            delt(i,j)=framR(i,j,1)-I_p(i,j);
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%% DCT of Delta matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    for i = 1:8:MM
        for j = 1:8:NN
           K = delt(i:i+7,j:j+7);
           K_DCT = dct2(K); 
           
           deltdct(i:i+7,j:j+7) = K_DCT;
        end
    end
    
    deldctfram(:,:,k)=deltdct(:,:);
end

%%%%%%%%%%%%%%%% Addition of DCT matrix of Delta and Frame %%%%%%%%%%%%%%%%

for i=1:6
    ffrmdct(:,:,i)=frmdct(:,:,i)-deldctfram(:,:,i);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
no_of_zeros=0;

for k=1:6;
    
    I_G=framR(:,:,k);
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
           K_DCT(abs(K_DCT)<(50)) = 0;
        
            for ii = 1:8
                for jj = 1:8
                    if(abs(K_DCT(ii,jj))<(50))
                        no_of_zeros=no_of_zeros+1;
                    end
                end
            end
        
            K_RE = idct2(K_DCT);
            I_re(i:i+7,j:j+7) = K_RE;
        
        end
    end
 refram(:,:,k)=I_re;    
 figure,
 imshow(uint8(refram(:,:,k)));

 figure,imshow(uint8(I_p));
 figure, imshow(uint8(refram(:,:,k))- uint8(I_p));

end
% err = (uint8(I_re) - I_p);
% err = err.*err;
% MSE(k+1) = sum(sum(err))/(MM*NN)
% Comp_ratio(k+1) = ((MM*NN) - no_of_zeros)/(MM*NN);
