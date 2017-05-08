%% Author     : Isidro C. Medina Jr.
%
% Description : A function that extracts relevent genes from a dataset 
%               given a list of genes. 
%
% Copyright   : (c) Isidro C. Medina Jr., March 2008
%% 
%  Inputs: data          data_1, data_2
%          genelist      list of genes & their RPV vlues 
%          top           number of genes to extract
%  Output: biomarkers    matrix of of potential biomarkers
%%
function [biomarkers] = relevant_genes(data_1, data_2, genelist, top)
 subset_1 = []; subset_2 = [];
  for m=1:top
	  x = genelist(m,1);                  % top genes based on RPV values
	  gene_1 = data_1(x,:);               % gene(m)samples from data
	  gene_2 = data_2(x,:);
	  subset_1 = cat(1,subset_1,gene_1);  % sample subset from data
	  subset_2 = cat(1,subset_2,gene_2);  
  end;
    subset = [subset_1,subset_2];       % combine expresion profiles from  
    idx = genelist(1:top,1);            % samples based on selected genes
    biomarkers = [idx, subset];          % pad biomarkers with gene index
end