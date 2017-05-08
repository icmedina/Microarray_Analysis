%% Author     : Isidro C. Medina Jr.
%
% Description : A function that determines the run times of several
%               processes
%
% Copyright   : (c) Isidro C. Medina Jr., March 2008
%%
% Input:  Run Times

%%
function [time] = run_times(time_var_1,time_var_2,time_p,time_n,time_np)
time1 = cat(1,time_var_1,time_var_2);
time2 = cat(1,time_p,time_n);
time3 = cat(1,time1,time2);
time4 = cat(1,time3,time_np);
sum_time = sum(time4)/60;
time = cat(1,time4,sum_time);
end