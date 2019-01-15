function [ E, dE ] = thetaEnergy( theta, q, p, d0, bCells, RBx, RBy, mask )
    % REDUCED ENERGY 
    
    dP = d0*p;
    dT = d0*theta;
    dQ = d0*q;
    QL = sqrt(sum(dQ.^2,2));
    
    rho = bsxfun(@rdivide,d0*bsxfun(@times,p,q),dP);
    Rsq = ( p(bCells(:,1)).*p(bCells(:,2)).*QL.^2 - (dP .* dT) ) ./ dP.^2;
    R = sqrt(Rsq);
    indZ = Rsq<0;
    R(indZ) = 0;
    
    deltaX = -bsxfun(@minus,RBx,rho(:,1));
    deltaY = -bsxfun(@minus,RBy,rho(:,2));
    dMag = sqrt(deltaX.^2 + deltaY.^2);
    d = bsxfun(@minus,dMag,R);

    E = .5*mean( sum(mask.*(d.^2),2) ); 
    if (nargout == 2)
        
        [ dR ] = fitDual.AFN.returnBondGradsTheta( q, theta, p, bCells );
        
%         dNormX = sum( (deltaX .* (d ./ dMag)) .* mask,2);
%         dNormY = sum( (deltaY .* (d ./ dMag)) .* mask,2);
        avgD = sum(d .* mask,2);
        NB = size(d0,1);
%         NC = size(q,1);
        II = 1:NB;
        indZ = ~indZ;
        dE = -sum( (sparse(II,II,avgD.*indZ)*dR) ,1)' / NB;
%         dE = dE((2*NC+1):3*NC);
    end
    
end

