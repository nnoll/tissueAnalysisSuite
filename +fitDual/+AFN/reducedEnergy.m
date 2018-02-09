function [ E, dE ] = reducedEnergy( x, d0, t1, t2, r1, r2, bCells )
    % REDUCED ENERGY 
    
    q = x(:,1:2);
    p = x(:,3);

    B = d0*bsxfun(@times,q,p);
    deltaP = d0*p;

    B1 = B - bsxfun(@times,deltaP,r1);
    B2 = B - bsxfun(@times,deltaP,r2);

    L1 = sqrt(sum(B1.^2,2));
    B1 = bsxfun(@rdivide,B1,L1);
    L2 = sqrt(sum(B2.^2,2));
    B2 = bsxfun(@rdivide,B2,L2);

    E = .5*mean( dot(B1,t1,2).^2 + dot(B2,t2,2).^2 );

    if (nargout == 2)
        NB = size(d0,1);
        NC = size(d0,2);
        [ drX1, drY1, drX2, drY2 ] = fitDual.AFN.reducedLocalGrads( q, p, bCells, r1 ,r2 );

        IP1 = dot(B1,t1,2);
        IP2 = dot(B2,t2,2);

%         II = 1:length(L1);
% 
%         T1x = sparse(II,II,IP1.*t1(:,1));
%         T1y = sparse(II,II,IP1.*t1(:,2));
%         T2x = sparse(II,II,IP2.*t2(:,1));
%         T2y = sparse(II,II,IP2.*t2(:,2));
        
%         dE = sum( (T1x * drX1) + (T1y * drY1) + ...
%                   (T2x * drX2) + (T2y * drY2), 1)' / NB;

        dE = bsxfun(@times,IP1.*t1(:,1),drX1) + bsxfun(@times,IP1.*t1(:,2),drY1) + ...
             bsxfun(@times,IP2.*t2(:,1),drX2) + bsxfun(@times,IP2.*t2(:,2),drY2);
        dE = dE / NB;
%         cols = [1:NB,1:NB,1:NB,1:NB,1:NB,1:NB;];
        rows = [bCells(:,1);bCells(:,1)+NC;bCells(:,1)+2*NC;bCells(:,2);bCells(:,2)+NC;bCells(:,2)+2*NC];
        vals = dE(:);
        dE = accumarray(rows,vals',[3*NC,1]);

    end
end

%         dBX1 = bsxfun(@rdivide,drX1,L1) - ...
%         bsxfun(@times,B1(:,1)./L1,bsxfun(@times,drX1,B1(:,1))+bsxfun(@times,drY1,B1(:,2)));
% 
%         dBY1 = bsxfun(@rdivide,drY1,L1) - ...
%         bsxfun(@times,B1(:,2)./L1,bsxfun(@times,drX1,B1(:,1))+bsxfun(@times,drY1,B1(:,2)));
% 
%         dBX2 = bsxfun(@rdivide,drX2,L2) - ...
%         bsxfun(@times,B2(:,1)./L2,bsxfun(@times,drX2,B2(:,1))+bsxfun(@times,drY2,B2(:,2)));
% 
%         dBY2 = bsxfun(@rdivide,drY2,L2) - ...
%         bsxfun(@times,B2(:,2)./L2,bsxfun(@times,drX2,B2(:,1))+bsxfun(@times,drY2,B2(:,2)));

%         dE = sum(bsxfun(@times,dBX1,IP1.*t1(:,1)) + bsxfun(@times,dBY1,IP1.*t1(:,2)) + ...
%                  bsxfun(@times,dBX2,IP2.*t2(:,1)) + bsxfun(@times,dBY2,IP2.*t2(:,2)),1)/Nb;
