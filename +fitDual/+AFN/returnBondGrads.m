function [ dRhoX, dRhoY, dR ] = returnBondGrads( q, theta, p, bCells )
    % RETURN BOND GRADS 

    p1 = p(bCells(:,1));
    p2 = p(bCells(:,2));
    q1x = q(bCells(:,1),1);
    q2x = q(bCells(:,2),1);
    q1y = q(bCells(:,1),2);
    q2y = q(bCells(:,2),2);
    t1 = theta(bCells(:,1));
    t2 = theta(bCells(:,2));
    
    dRhoX = fitDual.AFN.deriv.rhoXGrad(p1,p2,q1x,q2x);
    dRhoY = fitDual.AFN.deriv.rhoYGrad(p1,p2,q1y,q2y);
    dR = fitDual.AFN.deriv.radiusGrad(p1,p2,q1x,q2x,q1y,q2y,t1,t2);
    
%     NB = size(bCells,1);
%     NC = size(q,1);
%     rows = zeros(8*NB,1);
%     cols = zeros(8*NB,1);
%     valX = zeros(8*NB,1);
%     valY = zeros(8*NB,1);
%     valR = zeros(8*NB,1);
%     
%     for n = 1:8
%         II = (n-1)*NB + (1:NB);
%         rows(II) = 1:NB;
%         if (n <= 4)
%             cols(II) = bCells(:,1) + (n-1)*NC;
%         else
%             cols(II) = bCells(:,2) + (n-5)*NC;
%         end
%         valX(II) = dRhoX(:,n);
%         valY(II) = dRhoY(:,n);
%         valR(II) = dR(:,n);
%     end
%     
%     dRhoX = sparse(rows(valX~=0),cols(valX~=0),valX(valX~=0),NB,4*NC);
%     dRhoY = sparse(rows(valY~=0),cols(valY~=0),valY(valY~=0),NB,4*NC);
%     dR = sparse(rows(valR~=0),cols(valR~=0),valR(valR~=0),NB,4*NC);

end

