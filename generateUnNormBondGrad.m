syms q1x q1y t1 p1 q2x q2y t2 p2 rx ry

rX = ((p1.*q1x - p2.*q2x) - (p1-p2).*rx) ;
rY = ((p1.*q1y - p2.*q2y) - (p1-p2).*ry) ;

Dx = gradient(rX,[q1x,q1y,p1,q2x,q2y,p2]);
Dy = gradient(rY,[q1x,q1y,p1,q2x,q2y,p2]);

Hx = hessian(rX,[q1x,q1y,p1,q2x,q2y,p2]);
Hy = hessian(rY,[q1x,q1y,p1,q2x,q2y,p2]);

matlabFunction(Dx,'File','rXGrad.m');
matlabFunction(Dy,'File','rYGrad.m');
matlabFunction(Hx,'File','rXHess.m');
matlabFunction(Hy,'File','rYHess.m');