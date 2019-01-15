function [ E, dE ] = energy( Q, d0, bCells, RBx, RBy, mask )
    % FORCE BALANCE ENERGY 
    
    q = Q(:,1:2);
    theta = Q(:,3);
    p = Q(:,4);
    
    % Distance from each edge assuming finite presure everywhere
%     dP = d0*p;
%     dT = d0*theta;
%     dQ = d0*q;
%     QL = sqrt(sum(dQ.^2,2));
    
    dP = p(bCells(:,1)) - p(bCells(:,2));
    dT = theta(bCells(:,1)) - theta(bCells(:,2));
    dQ = q(bCells(:,1),:) - q(bCells(:,2),:);
    QL = sum(dQ.^2,2);
    
    rho = bsxfun(@rdivide,d0*bsxfun(@times,p,q),dP);
    Rsq = ( p(bCells(:,1)).*p(bCells(:,2)).*QL - (dP .* dT) ) ./ dP.^2;
    R = sqrt(Rsq);
    indZ = Rsq<=0;
    R(indZ) = 0;
    
    deltaX = -bsxfun(@minus,RBx,rho(:,1));
    deltaY = -bsxfun(@minus,RBy,rho(:,2));
    dMag = sqrt(deltaX.^2 + deltaY.^2);
    d = bsxfun(@minus,dMag,R);

    E = .5*mean( sum(mask.*(d.^2),2) ); 
    
    if (nargout == 2)
        NB = size(d0,1);
        NC = size(d0,2);
        [ dRhoX, dRhoY, dR ] = fitDual.AFN.returnBondGrads( q, theta, p, bCells );
        dNormX = sum( (deltaX .* (d ./ dMag)) .* mask,2);
        dNormY = sum( (deltaY .* (d ./ dMag)) .* mask,2);
        avgD = sum(d .* mask,2);
        indZ = ~indZ;
        dR(~indZ,:) = 0;
%         dE = sum( (sparse(II,II,dNormX)*dRhoX) + (sparse(II,II,dNormY)*dRhoY) - (sparse(II,II,avgD.*indZ)*dR) ,1)' / NB;
        dE = (bsxfun(@times,dNormX,dRhoX) + bsxfun(@times,dNormY,dRhoY) - bsxfun(@times,avgD,dR)) / NB;
        rows = [bCells(:,1);bCells(:,1)+NC;bCells(:,1)+2*NC;bCells(:,1)+3*NC;bCells(:,2);bCells(:,2)+NC;bCells(:,2)+2*NC;bCells(:,2)+3*NC];
        vals = dE(:)';
        dE = accumarray(rows,vals,[4*NC,1]);
    end
    
end

