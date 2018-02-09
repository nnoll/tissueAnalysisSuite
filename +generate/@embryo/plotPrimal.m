function [ ] = plotPrimal( this, color, mode )
    % PLOT PRIMAL 
    
    if (nargin == 1 || ~strcmp(mode,'2D'))
        if (nargin == 1)
            patch('Faces',this.Mesh.ConnectivityList,'Vertices',this.Mesh.Points,'FaceColor',[1,1,1]);
        else
            patch('Faces',this.Mesh.ConnectivityList,'Vertices',this.Mesh.Points,'FaceVertexCData',color,'EdgeColor',[0,0,0]);
            shading faceted
        end
        axis equal
        colorbar()
    else
        [ph,th] = cart2sph(this.Mesh.Points(:,1),this.Mesh.Points(:,2),this.Mesh.Points(:,3));
        ph(ph<0) = - pi - ph(ph<0);
        ph(ph>0) = pi - ph(ph>0);
        th = pi/2-th;
        Tri_phi = ph(this.Mesh.ConnectivityList); 
        Tri_phi = [Tri_phi,Tri_phi(:,1)];
        Tri_phi = diff(Tri_phi,1,2);
        goodTris = (sum(abs(Tri_phi)>pi/2,2)) < 1;
        patch('Faces',this.Mesh.ConnectivityList(goodTris,:),'Vertices',[ph,th],'FaceVertexCData',color,'EdgeColor',[0,0,0]);
        shading faceted
        axis([-pi,pi,0,pi])
        view([90,90])
    end

end

