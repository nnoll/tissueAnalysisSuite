function [X,Y] = computeCartesian( this, Z, P )  

    A = polyval(this.a_n,Z,this.aS,this.aMu);
    B = polyval(this.b_n,Z,this.bS,this.bMu);
    X0 = polyval(this.x_n,Z,this.xS,this.xMu);
    Y0 = polyval(this.y_n,Z,this.yS,this.yMu);
    
    X = A .* cos(P) + X0;
    Y = B .* sin(P) + Y0;

end