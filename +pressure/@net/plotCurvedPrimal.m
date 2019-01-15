function [ ] = plotCurvedPrimal( this, rV, color, mode )
    %PLOT PRIMAL 

    if (nargin >= 2 && ~isempty(rV))
        r = this.computePrimalVerts(rV);
    else
        r = this.computePrimalVerts();
    end
    
    if (nargin <= 2)
       color = 'b';
    end
    
    if (nargin <= 3)
        mode = 0;
    end
    
    [ rho, R, edgeVerts, ~, badEdges ] = this.computeEdgeArcs();

    if (nargin <= 2)
        mode = 0;
    end
    
    if ( mode == 1 )
        T = this.returnTension();
        Tmax = prctile(T,90);
        Tmin = prctile(T,2);
        T = (T-Tmin)/(Tmax-Tmin);
        T(T>1) = 1;
        T(T<0) = 0;
        
        cmap = hot(256);
        xR = linspace(0,1,256);
        Tcolor(:,1) = interp1(xR,cmap(:,1),T);
        Tcolor(:,2) = interp1(xR,cmap(:,2),T);
        Tcolor(:,3) = interp1(xR,cmap(:,3),T);
        good = 1:size(Tcolor,1);
        good(badEdges) = [];
    end
   
    hold on
    for b = 1:size(edgeVerts,1)
        if (isnan(R(b)) || isinf(R(b)))
            if (mode == 0)
                plot(r(edgeVerts(b,:),1),r(edgeVerts(b,:),2),'LineWidth',2,'Color',color)
            else
                plot(r(edgeVerts(b,:),1),r(edgeVerts(b,:),2),'LineWidth',2,'Color',Tcolor(good(b),:)) 
            end
        else
            theta = atan2(r(edgeVerts(b,:),2)-rho(b,2),r(edgeVerts(b,:),1)-rho(b,1));
            theta(theta<0) = theta(theta<0) + 2*pi;
            theta = sort(theta);
            if (theta(2) - theta(1) > pi)
                theta(2) = theta(2)-2*pi;
                theta = theta(2:-1:1);
            end
            theta = linspace(theta(1),theta(2),30);
  
            x = rho(b,1) + R(b)*cos(theta);
            y = rho(b,2) + R(b)*sin(theta);
            
            if (mode == 0)
                plot(x,y,'LineWidth',2,'Color',color)
            else
                plot(x,y,'LineWidth',2,'Color',Tcolor(good(b),:))
            end
            
        end
    end
    
    if (mode == 0)
        scatter(r(:,1),r(:,2),'b','filled')
    end
    
end

