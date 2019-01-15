function [ ] = plotTension( this )
    %PLOTDUAL 

    Xd = this.Mesh.circumcenter;
    q = this.Mesh.Points;
    T = this.d0 * q; 
    T = sqrt( sum(T.^2,2) );
    T = (T - min(T))/(max(T)-min(T));
    cmap = hot(256);
    x = linspace(0,1,256);
    Tcolor(:,1) = interp1(x,cmap(:,1),T);
    Tcolor(:,2) = interp1(x,cmap(:,2),T);
    Tcolor(:,3) = interp1(x,cmap(:,3),T);

    patch('Faces',this.Mesh.ConnectivityList,'Vertices',q,'FaceColor','w','EdgeColor','none')
    hold on
    for e = 1:size(this.d1,2)
        dualVerts = (this.d1(:,e) ~= 0);
        line(Xd(dualVerts,1),Xd(dualVerts,2),Xd(dualVerts,3),'Color',Tcolor(e,:),'LineWidth',2);
    end
%     scatter3(Xd(:,1),Xd(:,2),Xd(:,3),'k','filled')
    axis equal
end

