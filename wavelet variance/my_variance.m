%% Author     : Isidro C. Medina Jr.
%
% Description : A function that computes for the global maximum-overlap
%               estimator of the wavelet variance
%
% Copyright   : (c) Isidro C. Medina Jr., March 2008
%%
% Input:  data     raw data
%         wave     'db1', 'db2', 'db4', 'db8'
% Output: S        Matrix of raw data padded via circular extension (right to left)
%         A,D      Matrix of wavelet coefficients (A = Approximations, D = Details)
%         D2       Matrix of squared wavelet details
%         Xl       Matrix of local wavelet variances (at each level)
%         Xg       Vector of global wavelet variances of each genes
%%
function [Xl,Xg,time_var] = my_variance(data, wave)
%data = dlmread('sample.csv', ',');
%wave = 'db1';
[M,N] = size(data);             % no. of samples (M) & no. of genes (N)   
clear v; t0 = cputime;

[Nw,lvl,S] = my_resize(N,data); % resize data to become a multiple of 2^n

% SWT of genes across samples using maximum level of decomposition (lvl)
fprintf(1,'Calculating the maximum-overlap estimator of wavelet variance... ');
Xg = []; Xl = [];
for n=1:M                       % (at different levels)
X = [];                             
 for j=1:lvl
  % compute for lambda (l) and L lambda (t) at level j    
    [A,D] = swt(S(n,:),j,wave);
    D2 = D.*D;
    l = 2^j;                    
    [t] = my_Llambda(wave,l);
  
  % compute for the wavelet variance (vy2)
    vy2 = sum(D2(j,t:Nw))/(2*l*Nw); 
    X = [X,vy2];
 end;
   Xl = cat(1,Xl,X);
   Xave = mean(Xl(n,:));
   Xg = cat(1,Xg,Xave);
end;

  t1 = cputime;
   fprintf(1,'Done.\n');
   time_var = t1-t0;
  fprintf('Running Time: %6.3f sec\n\n', time_var);
end