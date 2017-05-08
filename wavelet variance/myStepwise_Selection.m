%% Author     : Isidro C. Medina Jr.
%
% Description : A program that creates a Support Vector Machines classifier 
%               based on top ranking genes. It also computes for N-fold 
%               cross-validation accuracy, mean and standard deviation.

% Note:         Genes are fed to the classifier one by one. 
%
% Copyright   : (c) Isidro C. Medina Jr., March 2008
%
%% Forward Selection & Support Vector Machines
% Input:        Training data, Test Data, Classifier
% Output:       Gene Prediction Accuracy, Mean, Standard Deviation
% Required  [1] Bioinformatics toolbox for Support Vector Machines 
%  Toolbox:     Training (svmtrain) & Classifying (svmclassify)
%           [2] Statistical Toolbox for Cross Validation(crossvalind)
%%
function [Accuracy,time_svm] = myStepwise_Selection(data,celltype,cell,RPV,top,replicate)

 groups = ismember(celltype,cell);

 fprintf(1,'Testing the Accuracy of Gene Predictors... \n\n');
 clear v; t0 = cputime;
  
 Acc = []; 
 for rep=1:replicate
   fprintf('Iteration: %g \n\n',rep);
 
% Randomly select training & test sets.  
     [train, test] = crossvalind('holdOut',celltype,0.4,'classes',{'serous','endomtrd','clear','mucinous'}); 
     cp = classperf(groups);
   
   res = [];  
   for gene_idx=2:top+1
     gene_set = data(2:end,2:gene_idx);

% Train the SVM classifier using a linear kernel
 fprintf('Training SVM with %g gene(s)... ',gene_idx-1);
 svmStruct = svmtrain(gene_set(train,:),groups(train)); fprintf(1,'Done.\n');
   
% Classify the test set.
 fprintf('Testing SVM performance using %g gene(s)... ',gene_idx-1);
  classes = svmclassify(svmStruct,gene_set(test,:)); fprintf(1,'Done.\n\n');
   
% Evaluate the performance of the classifier.   
    classperf(cp,classes,test); 
    accuracy = cp.CorrectRate*100;
    res = cat(1,res,accuracy);
  end;
  Acc = [Acc,res];
  
 end;
t1 = cputime;
time_svm = t1-t0;
   
% Show rank & gene index on the output
 index=[];
 for gene_idx=1:top
    idx = [gene_idx,RPV(gene_idx,1)]; 
    index = cat(1,index,idx);
 end;

% compute for the mean of the prediction accuracy
 Mean = [];
 for gene_idx=1:top
    mean_acc = mean(Acc(gene_idx,:)); 
    Mean = cat(1,Mean,mean_acc);
 end; 
 
% compute for the standard deviation of the prediction accuracy
 StdDev = [];
 for gene_idx=1:top
    std_dev = std(Acc(gene_idx,:)); 
    StdDev = cat(1,StdDev,std_dev);
 end; 

% combine results
 Acc1 = [index, Acc]; 
 Acc2 = [Acc1, Mean]; 
 Accuracy = [Acc2,StdDev]; 
 
 % Plot the results
  %svmStruct = svmtrain(gene_set(train,:),groups(train),'showplot',true);
  %title(sprintf('Kernel Function: %s',func2str(svmStruct.KernelFunction)),'interpreter','none'); 
  %svmclassify(svmStruct,gene_set(test,:),'showplot',true);
   
 fprintf('Finished.\nRunning Time: %6.3f sec\n\n', time_svm);
end