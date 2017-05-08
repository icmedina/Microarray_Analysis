%% Author     : Isidro C. Medina Jr.
%
% Description : A function that extracts top ranking genes (p/n_rpv)
%               from the raw data
%
% Copyright   : (c) Isidro C. Medina Jr., March 2008
%%
% Input:  data          Data processed with bioconductor using GCRMA (d1,d2)
%         top           Number of genes that will be extracted
%         RPV           List of genes that will be extracted (p_rpv, n_rpv)
% Output: biomarkers    Expression values of selected genes
%
%%
function [biomarkers,p_markers,n_markers] = my_top_genes(data_1,data_2,top,RPV,pRPV,nRPV)
%p_rpv = dlmread('pRPV.csv', ',');      % list of genes with positive RPV
%n_rpv = dlmread('nRPV.csv', ',');      % list of genes with negative RPV 
%d1 =  dlmread('sample1.csv', ',');
%d2 =  dlmread('sample2.csv', ',');
%top = 5;                               % number of genes that will be extracted
  clear v; t0 = cputime;

% Extract Biomarkers
    fprintf(1,'Extracting Biomarkers... ');
  pmarkers = transpose(my_relevant_genes(data_1,data_2,pRPV,top));
  nmarkers = transpose(my_relevant_genes(data_1,data_2,nRPV,top));
  markers = transpose(my_relevant_genes(data_1,data_2,RPV,top));
  
% Assign class to samples. data 1 = 0; data 2 = 1
 [m,n] = size (data_1);      
 id0 = [];  id0(n+1,1) = 0;

 [m,n] = size (data_2);
 id1 = [];
  for i=1:n
   id1 =cat(1,id1,1);
  end;    
 id = cat(1,id0,id1);
 
  p_markers = [id,pmarkers];
  n_markers = [id,nmarkers];
  biomarkers = [id,markers];  % combine markers from nRPV & pRPV
  
    t1 = cputime;
    fprintf(1,'Voila!\n');
    fprintf('%g Biomarkers were successfully extracted.\n', top);
    fprintf('Running Time:%6.3f sec\n\n', t1-t0);

  %dlmwrite('biomarkers.csv',biomarkers);
end
 