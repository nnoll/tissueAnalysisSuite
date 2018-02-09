N = 11;
alpha = .03;
beta = .8;
gamma = .01;
mode = 3;

[ q, ~, b0 ] = generate.meshGrid(N,alpha);
q(b0,:) = 1.05*q(b0,:);
p = ones(size(q,1),1);
z = beta*( rand(size(p)) - .5 ); 
p = p + z;
% p(b0) = 1;
theta = .2 * ( rand(size(p))) ;

[ tri, L, q, p, theta, contin ] = pressure.calculateTri( q, p, theta, b0 );
clear z b0

PN = pressure.net(q,theta,p,tri);
Struct = PN.returnStruct(gamma);
