function [ e_phi, e_z ] = computeSurfaceBasis( this, Z, P )
    
    A = polyval(this.a_n,Z,this.aS,this.aMu);
    B = polyval(this.b_n,Z,this.bS,this.bMu);
    DA = polyval(this.Da,Z,this.aS,this.aMu)./this.aMu(2);
    DB = polyval(this.Db,Z,this.bS,this.bMu)./this.bMu(2);
    DX = polyval(this.Dx,Z,this.xS,this.xMu)./this.xMu(2);
    DY = polyval(this.Dy,Z,this.yS,this.yMu)./this.yMu(2);
    
    X1 = -A.*sin(P);
    Y1 = B.*cos(P);
    R1 = sqrt(X1.^2 + Y1.^2);
    X1 = X1./R1; 
    Y1 = Y1./R1;
    
    e_phi = vertcat(X1(:)',Y1(:)');
    e_phi = vertcat(e_phi,zeros(1,size(e_phi,2)));
    
    X2 = DA .* cos(P) + DX;
    Y2 = DB .* sin(P) + DY;
    Z2 = ones(size(X2));
    R2 = sqrt(X2.^2 + Y2.^2 + 1);
    X2 = X2./R2;
    Y2 = Y2./R2;
    Z2 = Z2./R2;
    
    e_z = vertcat(X2(:)',Y2(:)',Z2(:)');

    e_phi = e_phi';
    e_z = e_z';
    
end

