function [ chi, chiC, chiCV, kappa, L, Struct ] = returnCompatibility( N, alpha, beta, scale )
%RETURNCOMPATIBILITY Summary of this function goes here
%   Detailed explanation goes here

    contin = 1;
    while (contin)
        [ q, ~, b0 ] = generate.meshGrid(N,alpha);
        q(b0,:) = 1.05*q(b0,:);
        p = ones(size(q,1),1);
        z = beta*( rand(size(p)) - .5 ); 
        p = p + z;
        theta = .01 * ( rand(size(p))) ;

        [ tri, L, q, p, theta, contin ] = pressure.calculateTri( q, p, theta, b0, scale );
        clear z b0
    end

    PN = pressure.net(q,theta,p,tri);
%     [ ~, R ] = PN.computeEdgeArcs(1);
    r = PN.computePrimalVerts();
    R = sqrt(sum( (PN.d1'*r).^2, 2) );
    T = PN.returnTension();
    deltaP = abs(PN.d0*PN.p);
    kappa =  (deltaP.*R)./T;
    kappa = kappa(kappa>0);
    
    Struct = PN.returnStruct(0);
    Struct = seg.threefold_cell(Struct);
    
    [ chi, chiC, chiCV ] = isogonal.buildControlDistribution( Struct, 0 );
        
    chi = iqr(chi{1});
    chiC = iqr(chiC{1});
    chiCV = iqr(chiCV{1});
    
end

