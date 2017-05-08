% LOAD DATA & DETERMINE ITS DIMENSION
data = dlmread('sample.csv', ',');	%read data
sz = size(data);                    %dimension of data
c = min(sz);                        %determine column size
%l = length(data) ;                 %determine length of row

% [1] APPLY 1-D DWT TO EACH COLUMN (SAMPLES)
    % Initialize Matrices
MA1 = [];                           %approx. coefs.
MD1 = [];                           %detail coefs.
MDT = [];                           %thresholded detail coefs.

    % 1-D DWT of a the genes using single-level decomposition (db1)
for n=1:c
    [cA1,cD1] = dwt(data(:,n),'db1'); %,'mode','sym'); cA = approximation coef; cD = detail coef
    MA1 = [MA1,cA1];                %construct matrix consisting of approxiamation coefficients
    MD1 = [MD1,cD1];                %construct matrix consisting of detail coefficients

% [2] APPLY THRESHOLDING TO HIGH FREQUENCY WAVELET COEFFIECIENTS (per column)
    % calculate the default threshold parameters using ddencmp
    [thr,sorh,keepapp] = ddencmp('den','wv',cD1);   %den=denoising; wv=wavepackets 
    MDD1 = wthresh(cD1,sorh,thr);                   %MDD1 = thresholded sample
    MDT =[MDT,MDD1];
end;

% [3] CREATE A MATRIX OF NON-ZERO HIGH FREQUENCY WAVELET COEFFIECIENTS 
x = length(MDT);                    % length of high freq coef
Idx=[];                             % Indices of rows where sum is not equal to zero
MDT1 = [];                          % Initialize matrix of non-zero values
for m=1:x
    total = sum(abs(MDT(m,:)));     % Determine the sum of the absolute values of each row
     if (total~=0)                  % if the sum is not equal to zero
        Idx =cat(1,Idx,m);          % populate the list of non-zero of indices
        N = MDT(m,:);           
        MDT1 = cat(1,MDT1,N);       % populate the matrix of nonzero values with non-zero row
     end;
end;

% [4] CONCATENATE MATRICES (low freq. & thr high freq.)
M = cat(1, MA1, MDT1);       

% [5] Export data to comma delimited files
dlmwrite('dwt.csv',M);              %save the combined matrices (comma delimited)
%dlmwrite('indx_db1.1.csv',Idx);    %save the indices (comma delimited)

%De-noise signal using global thresholding with threshold best basis. 
%MDDC1 = wdencmp('gbl',MD1,'db8',1,thr,sorh,keepapp);    % de-noising and compresssion

%Display the approximation and detail of a 1 level decomposition.
subplot(2,2,1); plot(data(2500:3000)); title('Original OS Data')
subplot(2,2,2); plot(MA1(2500:3000)); title('Approximations of OS 1-D DWT')
subplot(2,2,3); plot(MD1(2500:3000)); title('Details of OS 1-D DWT')
subplot(2,2,4); plot(MDT); title('Threshold Details OS 1-D DWT')