%% Author     : Isidro C. Medina Jr.
%
% Description : Adjust the number of samples for SWT if N ~= 2^n 
%
% Copyright   : (c) Isidro C. Medina Jr., March 2008
%% 
% Input:  data      data processed with bioconductor using GCRMA 
%         N         Original size of the data
% Output: Nw        New sapmle size
%         S         Padded data
%         j         Maximum level based on new size
%%

function [Nw,j,S] = resize(N,data)
   m = mod(log(N),log(2));            % determine if N is divisible by 2^n
 if (m~=0)
      fprintf(1,'Resizing data to 2^n... ');
    j = log(N)/log(2);
    rem = m/log(2);
    j = j - rem + 1;
    Nw = 2^(j);                       % compute for new sample size
    diff = Nw-N-1;
   for d=0:diff
     x = N-d;
     data = [data(:,x+d),data];       % pad data from right to left
   end;
     S = data;   
 else
 % Return the original values if the sample size (N) is a multiple of 2^n
    Nw = N;
    j = log(N)/log(2);
    S = data;   
 end;
  fprintf(1,'Done.\n');
end