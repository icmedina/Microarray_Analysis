%% [1] Load data and global variances
clear all;
 fprintf(1,'Loading data... ');
  data_1 = load('OvCa-test.csv'); data_2 = load('normal-test.csv');
  Xg_1 = load('global_spectra_OvCa_db2.csv'); Xg_2 = load('global_spectra_Norm_db2.csv');
 fprintf(1,'Done.\n'); 
  whos;
top = 100;

%% [2] Calculate for the Relative Percent Variation 
 %RPV = (Xg_1 - Xg_2)/Xg_1 * 100
 [RPV,pRPV, nRPV] = my_rpv(Xg_1,Xg_2);

%% [3] Sort the genes according to RPV values
 pRPV = sortrows(pRPV,-2);
 nRPV = sortrows(nRPV,-2);
 RPV = sortrows(RPV,-2); 

%% [4] Extract Top Genes from the original data  
 [biomarkers,p_markers,n_markers] = my_top_genes(data_1,data_2,top,RPV,pRPV,nRPV);

%% [5] Compare local wavelet variances of sample 1 and sample 2
 %wavemenu;

%% [6] Output results to a tab delimited file
fprintf(1,'Saving results... ');
  %dlmwrite('pRPV.csv',pRPV); dlmwrite('nRPV.csv',nRPV);dlmwrite('RPV.csv',RPV);
  dlmwrite('p_markers-100.csv',p_markers); dlmwrite('n_markers-100.csv',n_markers);
  dlmwrite('biomarkers-200.csv',biomarkers);
fprintf(1,'Done.\n');  