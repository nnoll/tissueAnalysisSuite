syms q1x q1y t1 p1 q2x q2y t2 p2

rhoX = (p1.*q1x - p2.*q2x) ./ (p1 - p2);
rhoY = (p1.*q1y - p2.*q2y) ./ (p1 - p2);

Hx = hessian(rhoX,[q1x,q1y,t1,p1,q2x,q2y,t2,p2]);
Hy = hessian(rhoY,[q1x,q1y,t1,p1,q2x,q2y,t2,p2]);

matlabFunction(Hx,'File','rhoXHess.m');
matlabFunction(Hy,'File','rhoYHess.m');
