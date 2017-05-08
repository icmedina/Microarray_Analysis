%% Author     : Isidro C. Medina Jr.
%
% Description : A function that computes for Relative Percent Variation 
%               between two classes
%
% Copyright   : (c) Isidro C. Medina Jr., March 2008
%%
% Input:  X1     Global wavelet variance of sample 1
%         Y1     Global wavelet variance of sample 2
% Output: pRPV   List of genes with positive RPV
%         nRPV   List of genes with negative RPV
%%
function [RPV,pRPV,nRPV] = my_rpv(X1,Y1)
%X1 = dlmread('global_spectral_1.csv', ',');Y1 = dlmread('global_spectral_2.csv', ',');
m = length(X1); n = length(Y1);

if (n~=m)
    disp('Warning number of genes in the samples did not match!');
else
      fprintf(1,'Calculating Relative Percent Variation... ');

Idx = [];           % gene index
pRPV = [];          % positive RPV
nRPV = [];          % negative RPV
RPV = [];           % combined RPV

 for i=1:m
  x = X1(i);
  y = Y1(i);
  rpv =(x - y)/x*100;   % compute for RPV
  Idx = cat(1,Idx,i);
    if (rpv > 0)        % identify whether rpv value is positive or not
      prpv = [Idx(i),rpv];   
      pRPV = cat(1,pRPV,prpv);
    else
      rpv =abs(y - x)/y*100;  
      nrpv = [Idx(i),rpv];     
      nRPV = cat(1,nRPV,nrpv);
    end;
 end;
 RPV = cat(1,pRPV,nRPV);
 fprintf(1,'Done.\n\n');

 %dlmwrite('pRPV.csv',pRPV);
 %dlmwrite('nRPV.csv',nRPV);
 end;
end