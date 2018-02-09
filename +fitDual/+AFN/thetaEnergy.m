function [ E, dE ] = thetaEnergy( theta, q, p, d0, bCells, RBx, RBy, mask )
    % REDUCED ENERGY 
    
%     dP = d0*p;
%     dT = d0*theta;
%     dQ = d0*q;
    dP = p(bCells(:,1))-p(bCells(:,2));
    dT = theta(bCells(:,1))-theta(bCells(:,2));
    dQ = q(bCells(:,1),:)-q(bCells(:,2),:);
    QL = sum(dQ.^2,2);
    
    rho = bsxfun(@rdivide,d0*bsxfun(@times,p,q),dP);
    Rsq = ( p(bCells(:,1)).*p(bCells(:,2)).*QL - (dP .* dT) ) ./ dP.^2;
%     Rsq([1015,1330,1537])
%     pause
    R = sqrt(Rsq);
    indZ = Rsq<=0;
    R(indZ) = 0;
    
    deltaX = -bsxfun(@minus,RBx,rho(:,1));
    deltaY = -bsxfun(@minus,RBy,rho(:,2));
    dMag = sqrt(deltaX.^2 + deltaY.^2);
    d = bsxfun(@minus,dMag,R);

    E = .5*mean( sum(mask.*(d.^2),2) );
    
    if (nargout == 2)
        
        [ dR ] = fitDual.AFN.returnBondGradsTheta( q, theta, p, bCells );
        
%         cdfplot(R(sum(isinf(dR),2)==2))
%         pause
%         dR(1015,:)
%         R(1015)
%         tmp = sum(mask.*(d.^2),2);
%         tmp(1015)
%         pause
        
        avgD = sum(d .* mask,2);
        NB = size(d0,1);
        NC = size(d0,2);
%         indZ = ~indZ;
        dR(indZ,:) = 0;

        dE = -bsxfun(@times,avgD,dR) / NB;
        rows = [bCells(:,1);bCells(:,2)];
        dE = accumarray(rows,dE(:),[NC,1]);
%         find(isinf(dE))
%         find(isnan(dE))
%         isreal(dE)
%         plot(dE)
%         pause
    end
    
end


