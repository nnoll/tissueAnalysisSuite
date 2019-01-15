function [ H ] = returnReducedHessian( q, p, bCells, r1, r2, T1X, T1Y, T2X, T2Y, norm  )
    % RETURN REDUCED HESSIAN 
    
    if (nargin == 9 || norm == 1)
%         NB = size(bCells,1);
%         NC = size(q,1);

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
%         row = zeros(21*NB,1);
%         col = zeros(21*NB,1);
%         val1X = zeros(21*NB,1);
%         val1Y = zeros(21*NB,1);
%         val2X = zeros(21*NB,1);
%         val2Y = zeros(21*NB,1);

%         n = 1;
%         for ii = 1:6
%             for jj = ii:6
% 
%                 II = ((n-1)*NB) + (1:NB);
% 
%                 if (ii <= 3)
%                     row(II) = bCells(:,1) + (ii-1)*NC;
%                 else
%                     row(II) = bCells(:,2) + (ii-4)*NC;
%                 end
% 
%                 if (jj <= 3)
%                     col(II) = bCells(:,1) + (jj-1)*NC;
%                 else
%                     col(II) = bCells(:,2) + (jj-4)*NC;
%                 end
% 
%                 val1X(II) = T1X.*H1x(:,ii,jj);
%                 val1Y(II) = T1Y.*H1y(:,ii,jj);
%                 val2X(II) = T2X.*H2x(:,ii,jj);
%                 val2Y(II) = T2Y.*H2y(:,ii,jj);
% 
%                 n = n + 1;
%             end
%         end
% 
%         rows = [row(val1X~=0);row(val1Y~=0);row(val2X~=0);row(val2Y~=0)];
%         cols = [col(val1X~=0);col(val1Y~=0);col(val2X~=0);col(val2Y~=0)];
%         vals = [val1X(val1X~=0);val1Y(val1Y~=0);val2X(val2X~=0);val2Y(val2Y~=0)];
%         h = sparse(rows,cols,vals,3*NC,3*NC);
    else
%         NB = size(bCells,1);
        NC = size(q,1);
        rows = [bCells(:,1),bCells(:,1)+NC,bCells(:,2),bCells(:,2)+NC];
        cols = [bCells(:,1)+2*NC,bCells(:,1)+2*NC,bCells(:,2)+2*NC,bCells(:,2)+2*NC];
        vals = [T1X,T1Y,-T2X,-T2Y];
        h = sparse(rows,cols,vals,3*NC,3*NC);

    end
end

