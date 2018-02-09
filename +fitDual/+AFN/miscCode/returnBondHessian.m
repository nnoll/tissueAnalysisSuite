function [ hRhoX, hRhoY, hR ] = returnBondHessian( q, theta, p, bCells )
    % RETURN BOND GRADS 
    
    %% Rho Hessian
    hRhoX = fitDual.AFN.deriv.rhoXHess(p(bCells(:,1)),p(bCells(:,2)),q(bCells(:,1),1),q(bCells(:,2),1));
    hRhoY = fitDual.AFN.deriv.rhoYHess(p(bCells(:,1)),p(bCells(:,2)),q(bCells(:,1),2),q(bCells(:,2),2));

    hRhoX = reshape(hRhoX,size(hRhoX,1),8,8);
    hRhoY = reshape(hRhoY,size(hRhoY,1),8,8);
    [ hRhoX ] = fitDual.AFN.deriv.loadHessian( real(hRhoX), bCells, size(q,1) );
    [ hRhoY ] = fitDual.AFN.deriv.loadHessian( real(hRhoY), bCells, size(q,1) );

    %% Radius Hessian.
    hR = fitDual.AFN.deriv.radiusHess(p(bCells(:,1)),p(bCells(:,2)),q(bCells(:,1),1),...
                 q(bCells(:,2),1),q(bCells(:,1),2),q(bCells(:,2),2),theta(bCells(:,1)),theta(bCells(:,2)));
    hR = reshape(hR,size(hR,1),8,8);
    [ hR ] = fitDual.AFN.deriv.loadHessian( real(hR), bCells, size(q,1) );
    
end

% %% Rho Hessian.
%     rows = zeros(7*NB,1);
%     colX = zeros(7*NB,1);
%     colY = zeros(7*NB,1);
%     valX = zeros(7*NB,1);
%     valY = zeros(7*NB,1);
%     
%     % q1, p1
%     n = 1;
%     II = ((n-1)*NB) + (1:NB);
%     rows(II) = 1:NB;
%     colX(II) = (bCells(:,1)) + (4*NC)*(3*NC + bCells(:,1) - 1);
%     colY(II) = (bCells(:,1)+NC) + (4*NC)*(3*NC + bCells(:,1) - 1);
%     valX(II) = -p(bCells(:,2))./dP.^2;
%     valY(II) = -p(bCells(:,2))./dP.^2;
%     
%     % q1, p2
%     n = 2;
%     II = ((n-1)*NB) + (1:NB);
%     rows(II) = 1:NB;
%     colX(II) = (bCells(:,1)) + (4*NC)*(3*NC + bCells(:,2) - 1);
%     colY(II) = (bCells(:,1)+NC) + (4*NC)*(3*NC + bCells(:,2) - 1);
%     valX(II) = -p(bCells(:,1))./dP.^2;
%     valY(II) = -p(bCells(:,1))./dP.^2;
% 
%     % q2, p1
%     n = 3;
%     II = ((n-1)*NB) + (1:NB);
%     rows(II) = 1:NB;
%     colX(II) = (bCells(:,2)) + (4*NC)*(3*NC + bCells(:,1) - 1);
%     colY(II) = (bCells(:,2)+NC) + (4*NC)*(3*NC + bCells(:,1) - 1);
%     valX(II) = -p(bCells(:,2))./dP.^2;
%     valY(II) = -p(bCells(:,2))./dP.^2;
% 
%     % q2, p2
%     n = 4;
%     II = ((n-1)*NB) + (1:NB);
%     rows(II) = 1:NB;
%     colX(II) = (bCells(:,2)) + (4*NC)*(3*NC + bCells(:,2) - 1);
%     colY(II) = (bCells(:,2)+NC) + (4*NC)*(3*NC + bCells(:,2) - 1);
%     valX(II) = -p(bCells(:,1))./dP.^2;
%     valY(II) = -p(bCells(:,1))./dP.^2;
% 
%     % p1, p1
%     n = 5;
%     II = ((n-1)*NB) + (1:NB);
%     rows(II) = 1:NB;
%     colX(II) = (bCells(:,1)+3*NC) + (4*NC)*(3*NC + bCells(:,1) - 1);
%     colY(II) = (bCells(:,1)+3*NC) + (4*NC)*(3*NC + bCells(:,1) - 1);
%     valX(II) = 2*(rho(:,1) - q(bCells(:,1),1))./dP.^2;
%     valY(II) = 2*(rho(:,2) - q(bCells(:,1),2))./dP.^2;
% 
%     % p1, p2
%     n = 6;
%     II = ((n-1)*NB) + (1:NB);
%     rows(II) = 1:NB;
%     colX(II) = (bCells(:,1)+3*NC) + (4*NC)*(3*NC + bCells(:,2) - 1);
%     colY(II) = (bCells(:,1)+3*NC) + (4*NC)*(3*NC + bCells(:,2) - 1);
%     valX(II) = ( q(bCells(:,1),1) + q(bCells(:,2),1) - 2*rho(:,1) )./dP.^2;
%     valY(II) = ( q(bCells(:,1),2) + q(bCells(:,2),2) - 2*rho(:,2) )./dP.^2;
%     
%     % q2, p2
%     n = 7;
%     II = ((n-1)*NB) + (1:NB);
%     rows(II) = 1:NB;
%     colX(II) = (bCells(:,2)+3*NC) + (4*NC)*(3*NC + bCells(:,2) - 1);
%     colY(II) = (bCells(:,2)+3*NC) + (4*NC)*(3*NC + bCells(:,2) - 1);
%     valX(II) = 2*(rho(:,1) - q(bCells(:,2),1))./dP.^2;
%     valY(II) = 2*(rho(:,2) - q(bCells(:,2),2))./dP.^2;
%     
%     hRhoX = sparse(rows,colX,valX,NB,(4*NC)^2);
%     hRhoY = sparse(rows,colY,valY,NB,(4*NC)^2);
