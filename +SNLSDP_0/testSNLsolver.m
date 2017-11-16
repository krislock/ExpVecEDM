%%*************************************************************************   
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
%% authors: Kim-Chuan Toh and Yinyu Ye
%%*************************************************************************
 
   close all;    

   eg = 3;
   nf = sqrt(10^(-20/10));  %% noise factor
   noisetype = 'multiplicative'; 'additive';  
   randstate = 1; 
%%
   rand('state',randstate); 
   if (eg==1)
      P0 = 0.30*[1,-1, 1, -1; 
                 1, 1, -1, -1];  %% anchor positions
      npts = 20; 
      dim  = 2; 
      PP   = rand(dim,npts)-0.5; %% sensor positions
      Radius   = 0.40;
      BoxScale = 100;  
      P0 = BoxScale*P0;  
      PP = BoxScale*PP; 
      Radius = BoxScale*Radius;
   elseif (eg==2)
      P0 = (1/3)*[1 2 2 1; 
                  1 2 1 2;
                  1 1 2 2]-0.5;  %% anchor positions
      npts = 50; 
      dim  = 3; 
      PP   = rand(dim,npts)-0.5; %% sensor positions
      Radius   = 0.5;
      BoxScale = 100;       
      P0 = BoxScale*P0;  
      PP = BoxScale*PP; 
      Radius = BoxScale*Radius;
   elseif (eg==3) 
      %% an example from molecular conformation
      P0 = []; 
      filename = 'pdb1GM2.txt';
      Porig = readPDB(filename); %% atom positions
      [dim,N] = size(Porig);
      center = Porig*ones(N,1)/N; 
      PP = Porig - center*ones(1,N); 
      BoxScale = 2*ceil(max(max(abs(PP))));
      Radius = 5;  
   end
   nfix = size(P0,2); 
   [dim,npts] = size(PP); 
%%
%% main 
%%
   fprintf('\n number of anchors = %2.0f',size(P0,2));
   fprintf('\n number of sensors = %2.0d',size(PP,2));
   fprintf('\n box scale         = %3.2f',BoxScale);
   fprintf('\n radius            = %3.2f',Radius);
   fprintf('\n %s noise, noise factor = %3.2e',noisetype,nf); 
%%
   OPTIONS.alpha       = 1; %% regularization parameter
   OPTIONS.refinemaxit = 1000; 
   OPTIONS.plotyes     = 1; 
   OPTIONS.PP          = PP;   
   OPTIONS.BoxScale    = BoxScale; 
   OPTIONS.nf          = nf; 
%%
   DD = randistance(P0,PP,Radius,nf,noisetype,randstate);
   [Xopt,Yopt] = SNLsolver(P0,DD,dim,OPTIONS);
   if (nfix==0); 
      tvar = max(0,diag(Yopt)'-sum(Xopt.*Xopt)); 
      Xopt = matchposition(Xopt,PP,tvar);
   end 
   errtrue = sum((Xopt-PP).*(Xopt-PP));   
   MSE     = sum(errtrue)/npts;
   MSEdB   = 10*log10(MSE); 
   fprintf('\n-------------------------------------------------------')
   fprintf('\n (noise factor)^2 = %3.1fdB, ',10*log10(nf^2)); 
   fprintf('\n mean square error (MSE) in estimated positions = %3.1fdB',MSEdB);
   fprintf('\n-------------------------------------------------------\n')
%%*************************************************************************
