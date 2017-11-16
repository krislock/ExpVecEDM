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
%% Input: 
%% P0: anchor position. If there is no anchor, input P0 = [];
%% DD: distance between pairs of points.
%%     DD = [sensor-sensor distance, anchor-sensor distance]. 
%% dim: dimensional of the sensor space.
%% OPTIONS: parameters (optional).
%%*************************************************************************

   function  [Xopt,Yopt,SDPInfo,RefineInfo] = SNLsolver(P0,DD,dim,OPTIONS)    

   refinemaxit = 1000;
   plotyes  = 0;  
   alpha    = []; 
   PP       = [];
   printyes = 1;
   if ~exist('OPTIONS'); OPTIONS = []; end

   if isfield(OPTIONS,'refinemaxit'); refinemaxit = OPTIONS.refinemaxit; end
   if isfield(OPTIONS,'plotyes'); plotyes = OPTIONS.plotyes; end
   if isfield(OPTIONS,'alpha');   alpha   = OPTIONS.alpha; end
   if (plotyes); PP = OPTIONS.PP; BoxScale = OPTIONS.BoxScale; end
   if isfield(OPTIONS,'printyes'); printyes = OPTIONS.printyes; end
   pars.gaptol = 1e-6; 
   pars.printlevel = 0;
   if isfield(OPTIONS,'gaptol');     pars.gaptol = OPTIONS.gaptol; end
   if isfield(OPTIONS,'printlevel'); pars.printlevel = OPTIONS.printlevel; end
%%   
%%   
   nfix = size(P0,2);
   npts = length(DD) - nfix;    
   if any(size(DD) ~= [npts,npts+nfix]); 
      error('dimension of DD is not compatiable with the number of sensors and anchors'); 
   end
   if (norm(tril(DD),'fro')==0); DD = triu(DD) + triu(DD)'; end;   
%%
%%
   tstart = clock;
   [blk,AA,cc,bb] = generateSDPdata(P0,DD,dim,alpha); 
   if (printyes)
      fprintf('\n-------------------------------------------------------')
      fprintf('\n        estimate sensor positions by SDP');
      fprintf('\n-------------------------------------------------------')
      fprintf('\n num of constraints = %4.0d,',length(bb));
      fprintf('\n Please wait:');
      fprintf('\n solving SDP by the SDPT3 software package\n');
   end
   [obj,xx,yy,zz,info] = sqlp(blk,AA,cc,bb,pars);
   sdpobj = bb'*yy;
   Zopt = xx{1};   
   SDPInfo = info; 
   Yopt = Zopt(1:npts,1:npts); 
   if (nfix==0)
      [V,Lam] = eig(full(Yopt)); [lam,ind] = sort(diag(Lam)); 
      lam = lam(end:-1:1); ind = ind(end:-1:1);
      Xopt = diag(sqrt(lam(1:dim)))*V(:,ind(1:dim))';
   else
      Xopt = Zopt(npts+[1:dim],1:npts);  
   end
   residual = zeros(1,npts); 
   for k = 1:npts
      ddtmp = DD([1:npts],k)'; 
      idx = find(ddtmp);
      tmp = Xopt(:,idx)-Xopt(:,k)*ones(1,length(idx)); 
      errtmp = sqrt(sum(tmp.*tmp)) - ddtmp(idx); 
      residual(1,k) = sqrt(sum(errtmp.^2))/sqrt(length(idx)); 
   end
   ttime = etime(clock,tstart); 
   if (info.termcode > 0); 
      fprintf('\n ***** SDP is infeasible *****\n'); 
   end
   if (printyes); fprintf(' sdpobj = %4.3e, time = %4.1fs\n',sdpobj,ttime); end
%%
%% plots
%%
   if (plotyes)
      tvar = max(0,diag(Yopt)'-sum(Xopt.*Xopt));   
      if (nfix==0); 
         Xtmp = matchposition(Xopt,PP,tvar);
      else
         Xtmp = Xopt;
      end    
      errtrue = sum((Xtmp-PP).*(Xtmp-PP));  
      RMSD = sqrt(sum(errtrue))/sqrt(npts); 
      fprintf(' RMSD = %4.2e',RMSD);
      figure(1)
      plotpositions(P0,PP,Xtmp,[],BoxScale);
      xlabel(['SDP: RMSD = ',sprintf('%4.2e',RMSD)]);
      title(['nf = ',num2str(OPTIONS.nf),',  \lambda = ',...
      sprintf('%3.1e',alpha)]); 
   end  
%%
   RefineInfo = [];
   if (refinemaxit)
      fprintf('\n-------------------------------------------------------')
      fprintf('\n        refine positions by steepest descent');
      fprintf('\n-------------------------------------------------------')
      tstart = clock; 
      plotidx = [];
      [XSD,Info] = refinepositions(Xopt,P0,DD,refinemaxit);
      Xopt = XSD;
      RefineInfo = Info; 
      ttime = etime(clock,tstart);
      len = length(Info.objective);
      if (printyes)
         fprintf('\n objstart = %5.4e, objend = %5.4e',...
         Info.objective(1),Info.objective(len));
         fprintf('\n number of iterations = %2.1d, time = %4.1fs\n',len,ttime);
      end     
      if (plotyes)
         if (nfix==0); Xopt = matchposition(Xopt,PP,tvar); end 
         errtrue = sum((Xopt-PP).*(Xopt-PP));   
	 RMSD = sqrt(sum(errtrue))/sqrt(npts); 
         fprintf(' RMSD = %4.2e',RMSD);
         figure(2)
         plotpositions(P0,PP,Xopt,[],BoxScale); 
         xlabel(['Refinement: RMSD = ',sprintf('%4.2e',RMSD)]);
         title(['nf = ',num2str(OPTIONS.nf),',  \lambda = ',...
         sprintf('%3.1e',alpha)]); 
      end
   end
%%*************************************************************************
