function [ H ] = returnReducedHessian( q, p, bCells, r1, r2, T1X, T1Y, T2X, T2Y, norm  )
    % RETURN REDUCED HESSIAN 
    
    if (nargin == 9 || norm == 1)
        
        p1 = p(bCells(:,1));
        p2 = p(bCells(:,2));
        q1x = q(bCells(:,1),1);
        q1y = q(bCells(:,1),2);
        q2x = q(bCells(:,2),1);
        q2y = q(bCells(:,2),2);

        H1x = fitDual.AFN.deriv.sXHess(p1,p2,q1x,q2x,q1y,q2y,r1(:,1),r1(:,2));
        H1y = fitDual.AFN.deriv.sYHess(p1,p2,q1x,q2x,q1y,q2y,r1(:,1),r1(:,2));
        H2x = fitDual.AFN.deriv.sXHess(p1,p2,q1x,q2x,q1y,q2y,r2(:,1),r2(:,2));
        H2y = fitDual.AFN.deriv.sYHess(p1,p2,q1x,q2x,q1y,q2y,r2(:,1),r2(:,2));

        H1x = bsxfun(@times,reshape(H1x,size(H1x,1),6,6),T1X);
        H1y = bsxfun(@times,reshape(H1y,size(H1y,1),6,6),T1Y);
        H2x = bsxfun(@times,reshape(H2x,size(H2x,1),6,6),T2X);
        H2y = bsxfun(@times,reshape(H2y,size(H2y,1),6,6),T2Y);

        H = H1x + H1y + H2x + H2y;

    else
        
        NB = size(bCells,1);
        H = zeros(NB,6,6);
        
        % Cell 1
        H(:,1,3) = T1X;
        H(:,3,1) = T1X;
        H(:,2,3) = T1Y;
        H(:,3,2) = T1Y;
        
        % Cell 2
        H(:,4,6) = -T2X;
        H(:,6,4) = -T2X;
        H(:,5,6) = -T2Y;
        H(:,6,5) = -T2Y;
        
    end
end

