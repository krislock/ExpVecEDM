%%************************************************************************
%% Plots anchor points, estimations of unknown points and actual locations
%% and the discrepancies between them
%% Blue Diamond : Anchor
%% Red Star : Estimated position of unknown point
%% Green Circle : Actual position of unknown point
%%
%% Input: 
%% P0: Anchor positions. If there is no anchor, input P0 = [];
%% PP: Actual positions of unknown points
%% X_opt: Estimated positions of unknown points
%% Plane 'xy','yz','xz' : Desired 2-D plane if points are 3-D
%%       'xyz'          : 3-D 
%%************************************************************************

   function plotpositions(P0,PP,Xopt,plane,BoxScale);

   if ~exist('BoxScale'); BoxScale = 1; end
   if ~exist('plane'); plane = []; end

   axes('FontSize',14,'FontWeight','bold');
   markersize = 8; 

   dim = size(Xopt,1);
   if (dim == 2) & isempty(plane); plane = 'xy'; end
   if (dim == 3) & isempty(plane); plane = 'xyz'; end
%% 
   if strcmp(plane,'xy') | strcmp(plane,'xz') | strcmp(plane,'yz')
      if strcmp(plane,'xy')
         idx1 = 1; idx2 = 2; 
      elseif strcmp(plane,'xz')
         idx1 = 1; idx2 = 3;
      elseif strcmp(plane,'yz')
         idx1 = 2; idx2 = 3;
      end
      plot(Xopt(idx1,:),Xopt(idx2,:),'*r','markersize',markersize); 
      hold on; grid on
      plot([Xopt(idx1,:); PP(idx1,:)],[Xopt(idx2,:); PP(idx2,:)],'b')
      plot(PP(idx1,:),PP(idx2,:),'og','markersize',markersize); 
      if ~isempty(P0)
         h = plot(P0(idx1,:),P0(idx2,:),'bd','markersize',markersize);    
         set(h,'linewidth',3,'color','b');     
      end
      axis('square'); axis(0.6*BoxScale*[-1,1,-1,1]);
      pause(0.1); hold off
   elseif strcmp(plane,'xyz')
      plot3(PP(1,:),PP(2,:),PP(3,:),'og','markersize',markersize); 
      hold on; grid on; 
      plot3(Xopt(1,:),Xopt(2,:),Xopt(3,:),'.r','markersize',markersize); 
      plot3([Xopt(1,:); PP(1,:)],[Xopt(2,:); PP(2,:)],[Xopt(3,:); PP(3,:)],'b');
      if ~isempty(P0)
         h = plot3(P0(1,:),P0(2,:),P0(3,:),'d','markersize',markersize); 
         set(h,'linewidth',3,'color','b');     
      end   
      axis('square'); axis(0.6*BoxScale*[-1,1,-1,1,-1,1])
      pause(0.1); hold off
   end
%%************************************************************************


