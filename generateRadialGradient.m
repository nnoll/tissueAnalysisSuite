syms q1x q1y t1 p1 q2x q2y t2 p2

R = sqrt( (p1.*p2.*( (q1x-q2x).^2 + (q1y - q2y).^2 ) - (p1-p2).*(t1-t2))./(p1-p2).^2 );
H = gradient(R,[q1x,q1y,t1,p1,q2x,q2y,t2,p2]);
matlabFunction(H,'File','radiusGrad.m');
