function [ rV ] = computeCellVerts( this )
%COMPUTECELLVERTS 

    Xd = this.Mesh.circumcenter;
    [rV(:,1),rV(:,2)] = cart2sph(Xd(:,1),Xd(:,2),Xd(:,3));
    rV(:,2) = pi/2 - rV(:,2);
    
end

