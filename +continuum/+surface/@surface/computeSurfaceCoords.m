function [ z, phi ] = computeSurfaceCoords( this, v )
%COMPUTESURFACECOORDS

    z = v(:,3);
    
    A = polyval(this.a_n,z,this.aS,this.aMu);
    B = polyval(this.b_n,z,this.bS,this.bMu);
    
    X0 = polyval(this.x_n,z,this.xS,this.xMu);
    Y0 = polyval(this.y_n,z,this.yS,this.yMu);
    
    phi = atan2((v(:,2)-Y0)./B,(v(:,1)-X0)./A);

end

