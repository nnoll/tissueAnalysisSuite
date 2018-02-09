function Hr = radiusHessTheta(p1,p2,q1x,q2x,q1y,q2y,t1,t2)
%RADIUSHESSTHETA
%    HR = RADIUSHESSTHETA(P1,P2,Q1X,Q2X,Q1Y,Q2Y,T1,T2)

%    This function was generated by the Symbolic Math Toolbox version 5.10.
%    06-Jul-2017 12:09:34

t4 = p1-p2;
t5 = q1x-q2x;
t6 = q1y-q2y;
t7 = 1.0./t4.^2;
t8 = t1-t2;
t9 = t4.*t8;
t10 = t5.^2;
t11 = t6.^2;
t12 = t10+t11;
t15 = p1.*p2.*t12;
t13 = t9-t15;
t16 = t7.*t13;
t14 = 1.0./(-t16).^(3.0./2.0);
t17 = t7.*t14.*(1.0./4.0);
Hr = [t7.*t14.*(-1.0./4.0),t17,t17,-t17];
