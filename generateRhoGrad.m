syms q1x q1y t1 p1 q2x q2y t2 p2

rhoX = (p1.*q1x - p2.*q2x) ./ (p1 - p2);
rhoY = (p1.*q1y - p2.*q2y) ./ (p1 - p2);

Dx = gradient(rhoX,[q1x,q1y,t1,p1,q2x,q2y,t2,p2]);
Dy = gradient(rhoY,[q1x,q1y,t1,p1,q2x,q2y,t2,p2]);

matlabFunction(Dx,'File','rhoXGrad.m');
matlabFunction(Dy,'File','rhoYGrad.m');