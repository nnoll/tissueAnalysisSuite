function [ dR ] = returnBondGradsTheta( q, theta, p, bCells )
    % RETURN BOND GRADS 

    p1 = p(bCells(:,1));
    p2 = p(bCells(:,2));
    q1x = q(bCells(:,1),1);
    q2x = q(bCells(:,2),1);
    q1y = q(bCells(:,1),2);
    q2y = q(bCells(:,2),2);
    t1 = theta(bCells(:,1));
    t2 = theta(bCells(:,2));
    
    dR = fitDual.AFN.deriv.radiusGradTheta(p1,p2,q1x,q2x,q1y,q2y,t1,t2);
    
end

