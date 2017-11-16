%%**************************************************************
%% Plot connectivity graph
%%
%%**************************************************************

  function plotgraph(P0,PP,DD)

  axes('FontSize',14,'FontWeight','bold');
  markersize = 8; 


  [dim,nfix] = size(P0);  
  [dim,npts] = size(PP); 

  Ds = DD(1:npts,1:npts); 
  Da = DD(1:npts,npts+(1:nfix)); 

  r1 = 1; r2 = 2; 
  h = plot(PP(r1,1),PP(r2,1),'or','markersize',6);    
  set(h,'linewidth',3);
  hold on;      
  for j = 2:npts           
     idx = find(Ds(1:j,j)); 
     if ~isempty(idx)
        len = length(idx); 
        Pj = PP(:,j)*ones(1,len); 
        plot([Pj(r1,:); PP(r1,idx)],[Pj(r2,:); PP(r2,idx)],'g');
     end
     h = plot(PP(r1,j),PP(r2,j),'or','markersize',6);    
     set(h,'linewidth',3);      
  end
  if ~isempty(P0)
     for j = 1:nfix         
        idx = find(Da(:,j)); 
        len = length(idx); 
        Pj = P0(:,j)*ones(1,len); 
        h = plot([Pj(r1,:); PP(r1,idx)],[Pj(r2,:); PP(r2,idx)],'b');
        set(h,'linewidth',2);    
        h = plot(P0(r1,j),P0(r2,j),'bd','markersize',markersize);   
        set(h,'linewidth',3);      
     end
  end
  grid on
  axis('square');
  hold off
%%**************************************************************
