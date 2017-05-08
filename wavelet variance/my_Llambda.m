%% Author     : Isidro C. Medina Jr.
%
% Description : A function that determines the value of L_lambda based on 
%               the type of wavelet filter and scale value (lambda).
%               L is the length of the wavelet filter.
%
% Copyright   : (c) Isidro C. Medina Jr., March 2008
%% 
% Input:  lambda    Scale Value
%         wave      wavelet filters:'db1', 'db2', 'db4', 'db8'
% Output: L lambda
%
%%
function [Llamb] = Llambda(wave,lambda)
  switch (wave)
    case 'db1'  % at db1, L = 2
        Llamb = ((lambda-1)*(2-1)) + 1;
    case 'db2'  % at db2, L = 4
        Llamb = ((lambda-1)*(4-1)) + 1;
    case 'db4'  % at db4, L = 8
         Llamb = ((lambda-1)*(8-1)) + 1;
    case 'db8'  % at db8, L = 16
         Llamb = ((lambda-1)*(16-1)) + 1;         
    otherwise
        error('Invalid wavelet filter');
  end
end