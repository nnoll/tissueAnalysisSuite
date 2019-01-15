syms q1x q1y t1 p1 q2x q2y t2 p2

rhoX = (p1.*q1x - p2.*q2x) ./ (p1 - p2);
rhoY = (p1.*q1y - p2.*q2y) ./ (p1 - p2);
R = sqrt( (p1.*p2.*( (q1x-q2x).^2 + (q1y - q2y).^2 ) - (p1-p2).*(t1-t2))./(p1-p2).^2 );

Dx = gradient(rhoX,[t1,t2]);
Dy = gradient(rhoY,[t1,t2]);
Dr = gradient(R,[t1,t2]);

matlabFunction(Dx,'File','rhoXGradTheta.m');
matlabFunction(Dy,'File','rhoYGradTheta.m');
matlabFunction(Dr,'File','radiusGradTheta.m');

Hx = hessian(rhoX,[t1,t2]);
Hy = hessian(rhoY,[t1,t2]);
Hr = hessian(R,[t1,t2]);

matlabFunction(Hx,'File','rhoXHessTheta.m');
matlabFunction(Hy,'File','rhoYHessTheta.m');
matlabFunction(Hr,'File','radiusHessTheta.m');
