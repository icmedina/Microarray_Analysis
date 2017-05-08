rma = dlmread('RMA.csv', ',');  %read data
sz = size(rma);
c = min(sz);                    %determine column size
l = length(rma) ;               %determine length of row

% Multilevel wavelet decomposition of a signal
MA3 = [];                               %initialize matrix consisting of approx. coef.
MD3 = [];
for n=1:c
    [C,L] = wavedec(rma(:,n),3,'db1');  %A level 3 decomposition of the signal (using db1 wavelet)
    cA3 = appcoef(C,L,'db1',3);         %Extract the level 3 approximation coefficients from C
    cD3 = detcoef(C,L,3);
    MA3= [MA3,cA3];                     %construct matrix consisting of approxiamation coefficients
    MD3 = [MD3,cD3];
end;

% Thresholding
% calculate the default threshold parameters using ddencmp
[thr,sorh,keepapp] = ddencmp('den','wv',rma); %den=denoising; wv=wavepackets 

%De-noise signal using global thresholding with threshold best basis. 
%MAD3 = wdencmp('gbl',MA3,'db1',1,thr,sorh,keepapp);
%MDD3 = wdencmp('gbl',MD3,'db1',1,thr,sorh,keepapp);

%Combine threshold matrices
M3 = cat(1,MAD3,MDD3);

%Display the approximation and detail of a 3 level decomposition.
subplot(2,2,1); plot(rma); title('Original Data')
subplot(2,2,2); plot(MA3); title('1-D DWT db1.3 Approximations')
subplot(2,2,3); plot(MAD3); title('1-D DWT db1.3 Threshold Approximations')
subplot(2,2,4); plot(MDD3); title('1-D DWT db1.3 Threshold Details')

dlmwrite('db1.3.csv',M3);    %save to comma delimited file