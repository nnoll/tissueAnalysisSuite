syms q1x q1y t1 p1 q2x q2y t2 p2 rx ry

sX = ((p1.*q1x - p2.*q2x) - (p1-p2).*rx) / ...
     sqrt( (p1.*q1x - p2.*q2x - (p1-p2).*rx).^2 + (p1.*q1y - p2.*q2y - (p1-p2).*ry).^2 );
sY = ((p1.*q1y - p2.*q2y) - (p1-p2).*ry) / ...
     sqrt( (p1.*q1x - p2.*q2x - (p1-p2).*rx).^2 + (p1.*q1y - p2.*q2y - (p1-p2).*ry).^2 );

Dx = gradient(sX,[q1x,q1y,p1,q2x,q2y,p2]);
Dy = gradient(sY,[q1x,q1y,p1,q2x,q2y,p2]);

Hx = hessian(sX,[q1x,q1y,p1,q2x,q2y,p2]);
Hy = hessian(sY,[q1x,q1y,p1,q2x,q2y,p2]);

matlabFunction(Dx,'File','sXGrad.m');
matlabFunction(Dy,'File','sYGrad.m');
matlabFunction(Hx,'File','sXHess.m');
matlabFunction(Hy,'File','sYHess.m');