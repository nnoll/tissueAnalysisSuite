function [] = plotSurface( this )

    [Z,phi] = meshgrid(linspace(this.z_B(1),this.z_B(2),100),linspace(0,2*pi,100));
    [X,Y] = this.computeCartesian( Z, phi );
    surf(X,Y,Z,ones(size(Z)))
%     [N] = this.computeNormal( Z, phi );
%     [e_phi,e_z] = this.computeSurfaceBasis( Z, phi );
%     hold all
%     quiver3(X(:),Y(:),Z(:),N(1,:)',N(2,:)',N(3,:)','r');
%     quiver3(X(:),Y(:),Z(:),e_phi(1,:)',e_phi(2,:)',e_phi(3,:)','b');
%     quiver3(X(:),Y(:),Z(:),e_z(1,:)',e_z(2,:)',e_z(3,:)','k');
    
end
