function [ ] = plotTensionMap( this, Tmax )
    % PLOT TENSION MAP 
    
    T = this.returnTension();
    if (nargin == 1)
        Tmax = prctile(T,95);
        Tmin = prctile(T,5);
    end
    
    T = (T-Tmin)/(Tmax-Tmin);
    T(T>1) = 1;
    T(T<0) = 0;
    
    % Convert tension to color.
    cmap = hot(256);
    x = linspace(0,1,256);
    Tcolor(:,1) = interp1(x,cmap(:,1),T);
    Tcolor(:,2) = interp1(x,cmap(:,2),T);
    Tcolor(:,3) = interp1(x,cmap(:,3),T);
    
    r = this.computePrimalVerts();
    
    hold on
    for e = 1:size(this.d1,2)
        verts = find(this.d1(:,e)~=0);
        if (length(verts) == 2)
            plot(r(verts,1),r(verts,2),'Color',Tcolor(e,:),'LineWidth',2)
        end
    end
end

