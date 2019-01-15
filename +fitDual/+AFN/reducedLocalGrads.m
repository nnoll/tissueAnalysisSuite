function [ drX1, drY1, drX2, drY2 ] = reducedLocalGrads( q, p, bCells, r1 ,r2, norm )
%REDUCEDLOCALGRADS Summary of this function goes here
%   Detailed explanation goes here

    p1 = p(bCells(:,1));
    p2 = p(bCells(:,2));
    q1x = q(bCells(:,1),1);
    q1y = q(bCells(:,1),2);
    q2x = q(bCells(:,2),1);
    q2y = q(bCells(:,2),2);
    
    if (nargin == 5 || norm == 1)
        drX1 = fitDual.AFN.deriv.sXGrad(p1,p2,q1x,q2x,q1y,q2y,r1(:,1),r1(:,2));
        drY1 = fitDual.AFN.deriv.sYGrad(p1,p2,q1x,q2x,q1y,q2y,r1(:,1),r1(:,2));
        drX2 = fitDual.AFN.deriv.sXGrad(p1,p2,q1x,q2x,q1y,q2y,r2(:,1),r2(:,2));
        drY2 = fitDual.AFN.deriv.sYGrad(p1,p2,q1x,q2x,q1y,q2y,r2(:,1),r2(:,2));
    else
        drX1 = fitDual.AFN.deriv.rXGrad(p1,p2,q1x,q2x,r1(:,1));
        drY1 = fitDual.AFN.deriv.rYGrad(p1,p2,q1y,q2y,r1(:,2));
        drX2 = fitDual.AFN.deriv.rXGrad(p1,p2,q1x,q2x,r2(:,1));
        drY2 = fitDual.AFN.deriv.rYGrad(p1,p2,q1y,q2y,r2(:,2));
    end
 
%     NB = size(bCells,1);
%     NC = size(q,1);
%     rows = zeros(6*NB,1);
%     cols = zeros(6*NB,1);
%     vals1X = zeros(6*NB,1);
%     vals1Y = zeros(6*NB,1);
%     vals2X = zeros(6*NB,1);
%     vals2Y = zeros(6*NB,1);

%     for n = 1:6
%         II = (n-1)*NB + (1:NB);
%         rows(II) = 1:NB;
%         if (n <= 3)
%            cols(II) = bCells(:,1) + (n-1)*NC;
%         else
%            cols(II) = bCells(:,2) + (n-4)*NC; 
%         end
%         
%         vals1X(II) = D1x(:,n);
%         vals1Y(II) = D1y(:,n);
%         vals2X(II) = D2x(:,n);
%         vals2Y(II) = D2y(:,n);
%     end
%     
%     drX1 = sparse(rows,cols,vals1X,NB,3*NC);
%     drY1 = sparse(rows,cols,vals1Y,NB,3*NC);
%     drX2 = sparse(rows,cols,vals2X,NB,3*NC);
%     drY2 = sparse(rows,cols,vals2Y,NB,3*NC);

end