%%
%% Estimate sensor positions based on semidefinite programming relaxation
%% of the non-convex optimization problem (1), followed by steepest descent
%% with backtracking line search on smooth unconstrained problem (2). 
%%
%% (1) min sum_{i,j} abs(gij) + sum_{k,j} abs(gkj) 
%%     s.t norm(xi-xj)^2-dij^2 = gij
%%         norm(ak-xj)^2-dkj^2 = gjk 
%%
%% (2) min sum_{i,j} (norm(xi-xj)-dij)^2 +  sum_{k,j} (norm(xj-ak)-djk)^2
%%
%% For details, see
%% [1] P. Biswas, T.-C. Liang, K.-C. Toh, T.-C. Wang, and Y. Ye,    
%%     Semidefinite programming approaches for sensor network localization 
%%     with noisy distance measurements, 
%%     IEEE Transactions on Automation Science and Engineering, 
%%     3 (2006), pp. 360--371.
%%

To learn how to use this Matlab software for sensor network localization, 
run the following m-file in the directory SNLSDP-0: 

>> testSNLsolver

The software needs the SDPT3 software package (included).
The user may have to compile (only need to be done once)
the mex-files needed in SDPT3 in the sub-directory
SDPT3-4.0-beta as follows: 

>> cd SDPT3-4.0
>> Installmex  
>> cd .. 

