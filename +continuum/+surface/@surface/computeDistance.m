function [ d ] = computeDistance( this, r )  

    z = r(3,:);
    
    a = polyval(this.a_n,z,this.aS,this.aMu);
    b = polyval(this.b_n,z,this.bS,this.bMu);
    
    x0 = polyval(this.x_n,z,this.xS,this.xMu);
    y0 = polyval(this.y_n,z,this.yS,this.yMu);
 
    phi = atan2((r(2,:)-y0)./b,(r(1,:)-x0)./a);

    d = sqrt( ((r(1,:)-a.*cos(phi)-x0).^2) + ((r(2,:)-b.*sin(phi)-y0).^2) );

    %Compute signs - to be used for projection back to surface.
    indX = (r(1,:).^2 + r(2,:).^2) > ( (a.*cos(phi)+x0).^2 + (b.*sin(phi)+y0).^2 );
    d( indX ) = -d( indX );

end
