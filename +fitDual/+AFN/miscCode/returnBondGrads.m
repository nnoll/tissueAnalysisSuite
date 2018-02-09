function [ dRhoX, dRhoY, dR ] = returnBondGrads( q, theta, p, bCells, rho, R, dQ, QL )
    % RETURN BOND GRADS 

    NB = size(bCells,1);
    NC = size(q,1);
    
    dP = p(bCells(:,1)) - p(bCells(:,2));
    dT = theta(bCells(:,1)) - theta(bCells(:,2));

    %% Build rho matrix.
    rows = zeros(4*NB,1);
    colX = zeros(4*NB,1);
    colY = zeros(4*NB,1);
    valX = zeros(4*NB,1);
    valY = zeros(4*NB,1);
    
    % q(bCells1,:)
    n = 1;
    II = ((n-1)*NB) + (1:NB);
    rows(II) = 1:NB;
    colX(II) = bCells(:,1);
    colY(II) = bCells(:,1) + NC;
    valX(II) = p(bCells(:,1))./dP;
    valY(II) = p(bCells(:,1))./dP;
    
    % q(bCells2,:)
    n = 2;
    II = ((n-1)*NB) + (1:NB);
    rows(II) = 1:NB;
    colX(II) = bCells(:,2);
    colY(II) = bCells(:,2) + NC;
    valX(II) = -p(bCells(:,2))./dP;
    valY(II) = -p(bCells(:,2))./dP;
    
    % p(bCells1)
    n = 3;
    II = ((n-1)*NB) + (1:NB);
    rows(II) = 1:NB;
    colX(II) = bCells(:,1) + 3*NC;
    colY(II) = bCells(:,1) + 3*NC;
    valX(II) = (q(bCells(:,1),1) - rho(:,1))./dP;
    valY(II) = (q(bCells(:,1),2) - rho(:,2))./dP;
    
    % p(bCells2)
    n = 4;
    II = ((n-1)*NB) + (1:NB);
    rows(II) = 1:NB;
    colX(II) = bCells(:,2) + 3*NC;
    colY(II) = bCells(:,2) + 3*NC;
    valX(II) = -(q(bCells(:,2),1) - rho(:,1))./dP;
    valY(II) = -(q(bCells(:,2),2) - rho(:,2))./dP;

    dRhoX = sparse(rows,colX,valX,NB,4*NC);
    dRhoY = sparse(rows,colY,valY,NB,4*NC);

    %% Build radius matrix.
    rows = zeros(8*NB,1);
    cols = zeros(8*NB,1);
    vals = zeros(8*NB,1);
    
    % cell 1, qx
    n = 1;
    II = ((n-1)*NB) + (1:NB);
    rows(II) = 1:NB;
    cols(II) = bCells(:,1);
    vals(II) = (p(bCells(:,1)).*p(bCells(:,2)) .* dQ(:,1)) ./ (R.*dP.^2);
    
    % cell 1, qy
    n = 2;
    II = ((n-1)*NB) + (1:NB);
    rows(II) = 1:NB;
    cols(II) = bCells(:,1) + NC;
    vals(II) = (p(bCells(:,1)).*p(bCells(:,2)) .* dQ(:,2)) ./ (R.*dP.^2);
    
    % cell 2, qx
    n = 3;
    II = ((n-1)*NB) + (1:NB);
    rows(II) = 1:NB;
    cols(II) = bCells(:,2);
    vals(II) = -(p(bCells(:,1)).*p(bCells(:,2)) .* dQ(:,1)) ./ (R.*dP.^2);
    
    % cell 2, qy
    n = 4;
    II = ((n-1)*NB) + (1:NB);
    rows(II) = 1:NB;
    cols(II) = bCells(:,2) + NC;
    vals(II) = -(p(bCells(:,1)).*p(bCells(:,2)) .* dQ(:,2)) ./ (R.*dP.^2);
    
    % cell 1, Theta
    n = 5;
    II = ((n-1)*NB) + (1:NB);
    rows(II) = 1:NB;
    cols(II) = bCells(:,1) + 2*NC;
    vals(II) = -1 ./ (2*R.*dP);
    
    % cell 2, Theta
    n = 6;
    II = ((n-1)*NB) + (1:NB);
    rows(II) = 1:NB;
    cols(II) = bCells(:,2) + 2*NC;
    vals(II) = 1 ./ (2*R.*dP);
    
    % cell 1, p
    n = 7;
    II = ((n-1)*NB) + (1:NB);
    rows(II) = 1:NB;
    cols(II) = bCells(:,1) + 3*NC;
    vals(II) =  p(bCells(:,2)).*QL.^2./ (2*R.*dP.^2) - p(bCells(:,1)).*p(bCells(:,2)).*QL.^2./(R.*dP.^3) + dT./(2*dP.^2.*R);
    
    % cell 1, p
    n = 8;
    II = ((n-1)*NB) + (1:NB);
    rows(II) = 1:NB;
    cols(II) = bCells(:,2) + 3*NC;
    vals(II) = p(bCells(:,1)).*QL.^2./ (2*R.*dP.^2) + p(bCells(:,1)).*p(bCells(:,2)).*QL.^2./(R.*dP.^3) - dT./(2*dP.^2.*R);
    
    dR = sparse(rows,cols,vals,NB,4*NC);
    
end

