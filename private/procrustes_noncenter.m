function [rmsd, Pout] = procrustes_noncenter(Porig, P)
% [rmsd, Pout] = procrustes_noncenter(Porig, P)
%                    rmsd=(n^(-0.5))*norm(Porig-P,'fro');
%
% This function performs a Procrustes rotation on P to match with Porig 
% in the best possible way.
% 
% INPUT: 
% Porig : n-by-r matrix, giving the original location of points
% P     : n-by-r matrix, giving the estimated location of points
% **Note: Porig and P do NOT need to be centered before being input
% 
% OUPTUT:
% rmsd      : rmsd between Porig and P 
%             rmsd=n^(-0.5)*norm(Porig-P,'fro')
% Pout      : P after Procrustes rotation

% Updated: Sep 24, 2014 by Nathan
% Updated: May 13, 2016 by Henry
% Updated: Jun 20, 2016 by Nathan


n=size(P,1);

Porigshift=sum(Porig)/n;
Pshift=sum(P)/n;

[uu,ss,vv]=svd(( P-ones(n,1)*Pshift )'*( Porig-ones(n,1)*Porigshift ));
 % study the singular values and see if they are close to 1
 if min(diag(ss))<0.01,
     fprintf('min(ss)<0.01 in procrustes_noncenter \n');
     %keyboard
 end
optshift=Porigshift'-vv*(uu'*Pshift');

%Porig=Porig+ones(n,1)*Porigshift;
Pout=(P*uu)*vv'+ones(n,1)*optshift';

% calculate rmsd
rmsd=(n^(-0.5))*norm(Porig-Pout,'fro');

% % calculate "new" rmsd
% rmsd_new=0;
% for jj=1:n,
%     rmsd_new=rmsd_new+...
%         norm(Porig(jj,:)-P(jj,:))/(n*norm(Porig(jj,:)-Porigshift));
% end
