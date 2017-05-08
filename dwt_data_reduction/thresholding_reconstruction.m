%Construct approximations and details from the coefficients. 
%To construct the level 1 approximation and detail (A1 and D1) from the coefficients cA1 and cD1, type 
A1 = upcoef('a',cA,'db1',1,l); 
D1 = upcoef('d',cD,'db1',1,l);

%or using inverse wavelet transform,idwt
A1 = idwt(cA,[],'db1',l); 
D1 = idwt([],cD,'db1',l);

% Extract the levels 3, 2, and 1 detail coefficients from C
cD3 = detcoef(C,L,3);
cD2 = detcoef(C,L,2);
cD1 = detcoef(C,L,1);
% or 
[cD1,cD2,cD3] = detcoef(C,L,[1,2,3]);



%find default threshold values per column
th=0;
for n=1:c
    [thr,sorh,keepapp] = ddencmp('den','wv',MA(:,n)); %den=denoising; wv=wavepackets 
        if  (th<thr)
           th=thr;
        end;
end;