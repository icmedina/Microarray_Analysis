% LOAD DATA
n_rpv = dlmread('nRPV.csv', ',');
p_rpv = dlmread('pRPV.csv', ',');
s1 =  dlmread('sample1.csv', ',');
s2 =  dlmread('sample2.csv', ',');

% [1] Read list of coef
l_coef = length(top100_cf);             % nos. of genes (length)
M_gs = [];                              % initialize gene subset matrix
for m=1:l_coef
    coef1 = top100_cf(m);
        gs1 = data(coef1,:);
    coef2 = coef1 - 1;
        gs2 = data(coef2,:);
    M_gs = cat(1,M_gs,gs1);              % [2] populate gene subset matrix
    M_gs = cat(1,M_gs,gs2);
end;

% [3] create list of genes
top_genes = [];
for m=1:l_coef
    coef1 = top100_cf(m);
    coef2 = coef1 - 1;
    top_genes = cat(1,top_genes,coef1);
    top_genes = cat(1,top_genes,coef2);
end;

subset = [top_genes,M_gs];
dlmwrite('gene_subset.csv',subset);