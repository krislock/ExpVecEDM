%%*************************************************************************
%% Generate SDP constraints and objective for SDPT3 corresponding to 
%% Sum of absolute errors
%%   min sum_{k,j} abs(gkj) + sum_{i,j} abs(gij)
%%   s.t norm(ak-xj)^2-dkj^2 = gjk 
%%       norm(xi-xj)^2-dij^2 = gij
%%
%% Input: P0 = anchor positions
%%        DD = npts x (nfix+npts) 
%%           = [unknown point-unknown point distances, 
%%              anchor-unknown point distances]
%%        dim = dimension of space
%%        alpha : Regularization parameter
%%       
%% Output: blk: Number and type of variables
%%         At,bb: Constraint equation coefficients 
%%         C: Objective function coefficient
%%*************************************************************************

   function  [blk,At,C,bb] = generateSDPdata(P0,DD,dim,alpha); 

   if ~exist('alpha'); alpha = []; end

   nfix = size(P0,2);         %% Number of anchors
   npts = size(DD,2) - nfix;  %% Number of unknown senor points

   if any(size(DD) ~= [npts,npts+nfix]); 
      error('dimension of DD is not compatiable with the number of sensors and anchors'); 
   end
%%
   if isempty(alpha) 
     %% If no anchors are provided, set alpha = 1 to stretch points apart
     alpha = 1.0; 
   end 
   if (nfix==0); dim = 0; end  
%%
%% Fix lower 2x2 block to be identity matrix
%% Note: Z = [Y X'; X I]
%%
   mm = 0; %% Counter variable with number of constraints
   for i = 1:dim
      for j = i:dim
         mm = mm+1;
         zz = sparse(npts+dim,1); 
	 if (j==i)
            zz(npts+i,mm) = 1;  
            bb(mm,1) = 1;
         else
            zz(npts+[i,j]) = [1;1];
            bb(mm,1) = 2;
         end
         AA{1,mm} = zz*zz'; 
      end
   end  
   if (nfix == 0) 
      %% If no anchors are present, 
      %% set e'*Y*e = 0 to center points around origin
      mm = mm+1; 
      zz = [ones(npts,1); zeros(dim,1)];
      AA{1,mm} = zz*zz';
      bb(mm,1) = 0; 
   end
%%     
%% Set constraints corresponding to distance measures 
%%
   for k = 1:npts
      for j = k+1:npts+nfix
         rr = DD(k,j);
         if (rr)  % If distance measure exists 
            mm = mm+1;  
            zz = sparse(npts+dim,1); 
            if (j > npts)
                zz([k,npts+(1:dim)]) = [-1; P0(:,j-npts)]; % Anchor-unknown constraint
            else
                zz([k,j]) = [1;-1];  % Unknown-unknown constraint
            end
            AA{1,mm} = zz*zz'; 
            bb(mm,1) = rr^2;
         end
      end
   end
%%
%% Semidefinite constraints
%%
   blk{1,1} = 's'; blk{1,2} = npts+dim; 
   At = svec(blk,AA,1); 
%%
%% Nonegative Slack variables
%%
   blk{2,1} = 'l'; blk{2,2} = 2*mm;
   At{2,1} = [-speye(mm); speye(mm)]; 
   C{2,1} = ones(2*mm,1); 
%%
%% Regularization term
%%
   if (nfix > 0)
      e = [ones(npts,1); P0*ones(nfix,1)]/sqrt(npts+nfix); 
   else
      e = ones(npts,1)/sqrt(npts); 
   end
   C{1} = -alpha*(speye(npts+dim,npts+dim)-e*e');
%%************************************************************************
