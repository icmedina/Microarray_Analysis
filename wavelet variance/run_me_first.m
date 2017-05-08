%% Author     : Isidro C. Medina Jr.
%
% Description : A program that computes for the maximum-overlap estimator 
%               of two samples then calculates the relative percent 
%               variation (RPV) with respect to the first sample.
%
% Copyright   : (c) Isidro C. Medina Jr., March 2008
%% 
% Input:  data     data processed with bioconductor using GCRMA 
%         wave     wavelet filters:'db1', 'db2', 'db4', 'db8'
% Output: Xl       Matrix of local wavelet variances (at each level)
%         Xg       Vector of global wavelet variances of each genes
%%
clear all;
 fprintf(1,'Loading data... ');
 load ('data/OvCa.mat');
 cell = 'serous';
  label = ismember(celltype,cell);
  data_1 = eset(:,label); [m n] = size(data_1);
  data_2 = eset(:,n+1:end);
 fprintf(1,'Done.\n'); whos
 top = 100;     % number of genes that will be extracted
 trial = 10;     % number of cross validation fold
 wave = 'db1';
   
%% [1] Computes for the wavelet variance
 [Xl_1, Xg_1,time_var_1] = my_variance(data_1, wave);
 [Xl_2, Xg_2,time_var_2] = my_variance(data_2, wave);

%% [2] Calculate for the Relative Percent Variation 
 %RPV = (Xg_1 - Xg_2)/Xg_1 * 100
 [RPV,pRPV, nRPV] = my_rpv(Xg_1,Xg_2);

%% [3] Sort the genes according to RPV values
 pRPV = sortrows(pRPV,-2);
 nRPV = sortrows(nRPV,-2);
 RPV = sortrows(RPV,-2); 

%% [4] Extract Top Genes from the original data  
 [biomarkers,p_markers,n_markers] = my_top_genes(data_1,data_2,top,RPV,pRPV,nRPV);
 
%% [5] Caculate Classification Accuracy using SVM
[Accuracy_P,time_p] = myStepwise_Selection(p_markers,celltype,cell,pRPV,top,trial);
[Accuracy_N,time_n] = myStepwise_Selection(n_markers,celltype,cell,nRPV,top,trial);
[Accuracy,time_np] = myStepwise_Selection(biomarkers,celltype,cell,RPV,top,trial);

[time] = run_times(time_var_1,time_var_2,time_p,time_n,time_np); 

%% [6] Output results to a tab delimited file

fprintf(1,'Saving results... ');
  dlmwrite('lcl_var_1.csv',Xl_1);  dlmwrite('lcl_var_2.csv',Xl_2);
  dlmwrite('glbl_var_1.csv',Xg_1); dlmwrite('glbl_var_2.csv',Xg_2);
  dlmwrite('pRPV.csv',pRPV); dlmwrite('nRPV.csv',nRPV);dlmwrite('RPV.csv',RPV);
  dlmwrite('p_markers.csv',p_markers); dlmwrite('n_markers.csv',n_markers);
  dlmwrite('biomarkers.csv',biomarkers); dlmwrite('Accuracy.csv',Accuracy);  
  dlmwrite('Accuracy-P.csv',Accuracy_P); dlmwrite('Accuracy-N.csv',Accuracy_N);
  dlmwrite('RunTimes.csv',time);
fprintf(1,'Done.\n');
fprintf('Over all Running Time: %6.3f mins\n\n', time(6));

%% [7] Compare local wavelet variances of sample 1 and sample 2
 %wavemenu;
