function [ ] = edgeHM( Struct, T )
%SKELETONIZE_VERTS Takes the network topology and triangulated position of
%vertices to produce a skeletonized image.

    Vlist = Struct.Vdat;
    Tmax = prctile(T,98);
    % Tmax = max(T);
    Tmin = prctile(T,2);
    % Tmin = min(T);

    T = (T-Tmin)/(Tmax-Tmin);
    T(T>1) = 1;
    T(T<0) = 0;
    
    % Convert tension to color.
    cmap = hot(256);
    x = linspace(0,1,256);
    
    Tcolor(:,1) = interp1(x,cmap(:,1),T);
    Tcolor(:,2) = interp1(x,cmap(:,2),T);
    Tcolor(:,3) = interp1(x,cmap(:,3),T);
    
    [ ~, dV, ~, iverts ] = fitDual.ATN.computeDiffOperators( Struct, 1 );  
    hold on
    for b = 1:size(dV,1)
        v = iverts(dV(b,:)==1);
        nv = iverts(dV(b,:)==-1);
        plot([Vlist(v).vertxcoord,Vlist(nv).vertxcoord],[Vlist(v).vertycoord,Vlist(nv).vertycoord],'Color',Tcolor(b,:),'LineWidth',2);
    end
    
end

