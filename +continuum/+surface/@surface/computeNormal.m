function [ N ] = computeNormal( this, Z, P )  

    if (nargin == 2)
        rv = Z;
        Z = rv(3,:);
        A = polyval(this.a_n,Z,this.aS,this.aMu);
        B = polyval(this.b_n,Z,this.bS,this.bMu);
        X0 = polyval(this.x_n,Z,this.xS,this.xMu);
        Y0 = polyval(this.y_n,Z,this.yS,this.yMu);
        P = atan2((rv(2,:)-Y0)./B,(rv(1,:)-X0)./A);
    else
        A = polyval(this.a_n,Z,this.aS,this.aMu);
        B = polyval(this.b_n,Z,this.bS,this.bMu);
        X0 = polyval(this.x_n,Z,this.xS,this.xMu);
        Y0 = polyval(this.y_n,Z,this.yS,this.yMu);
    end

    DA = polyval(this.Da,Z,this.aS,this.aMu)./this.aMu(2);
    DB = polyval(this.Db,Z,this.bS,this.bMu)./this.bMu(2);
    DX = polyval(this.Dx,Z,this.xS,this.xMu)./this.xMu(2);
    DY = polyval(this.Dy,Z,this.yS,this.yMu)./this.yMu(2);
    
    X = B .* cos(P);
    Y = A .* sin(P);
    Z = -( (A.*DB.*(sin(P).^2) + A.*DY.*sin(P)) + (B.*DA.*(cos(P).^2) + B.*DX.*cos(P)) );

    R = sqrt(X.^2 + Y.^2 + Z.^2);
    X = X./R;
    Y = Y./R;
    Z = Z./R;

    N = horzcat(X(:),Y(:),Z(:))';
end

