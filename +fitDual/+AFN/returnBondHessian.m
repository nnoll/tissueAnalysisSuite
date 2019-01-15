function [ h ] = returnBondHessian( q, theta, p, bCells, dNormX, dNormY, dAvg, lambda )
    
    %% Calculculate elements of Rho and Radius Hessian
    hRhoX = fitDual.AFN.deriv.rhoXHess(p(bCells(:,1)),p(bCells(:,2)),q(bCells(:,1),1),q(bCells(:,2),1));
    hRhoY = fitDual.AFN.deriv.rhoYHess(p(bCells(:,1)),p(bCells(:,2)),q(bCells(:,1),2),q(bCells(:,2),2));
    hR = fitDual.AFN.deriv.radiusHess(p(bCells(:,1)),p(bCells(:,2)),q(bCells(:,1),1),...
                 q(bCells(:,2),1),q(bCells(:,1),2),q(bCells(:,2),2),theta(bCells(:,1)),theta(bCells(:,2)));
    hCon = fitDual.AFN.deriv.constHess( p(bCells(:,1)), p(bCells(:,2)), q(bCells(:,1),1),...
                                q(bCells(:,2),1), q(bCells(:,1),2), q(bCells(:,2),2));

    hRhoX = reshape(hRhoX,[size(hRhoX,1),8,8]);
    hRhoY = reshape(hRhoY,[size(hRhoY,1),8,8]);
    hR = reshape(hR,[size(hR,1),8,8]);
    hCon = reshape(hCon,[size(hCon,1),8,8]);

    h = bsxfun(@times,hRhoX,dNormX) + bsxfun(@times,hRhoY,dNormY) - ...
        bsxfun(@times,hR,dAvg) + bsxfun(@times,hCon,lambda);
    
%     %% Add together in sparse matrix.
%     NB = size(bCells,1);
%     NC = size(q,1);
%     
%     row = zeros(36*NB,1);
%     col = zeros(36*NB,1);
%     valX = zeros(36*NB,1);
%     valY = zeros(36*NB,1);
%     valR = zeros(36*NB,1);
%     valC = zeros(36*NB,1);
%     
%     n = 1;
%     for ii = 1:8
%         for jj = ii:8
%             
%             II = ((n-1)*NB) + (1:NB);
%             
%             if (ii <= 4)
%                 row(II) = bCells(:,1) + (ii-1)*NC;
%             else
%                 row(II) = bCells(:,2) + (ii-5)*NC;
%             end
%             
%             if (jj <= 4)
%                 col(II) = bCells(:,1) + (jj-1)*NC;
%             else
%                 col(II) = bCells(:,2) + (jj-5)*NC;
%             end
%             
%             valX(II) = dNormX.*hRhoX(:,ii,jj);
%             valY(II) = dNormY.*hRhoY(:,ii,jj);
%             valR(II) = -dAvg.*hR(:,ii,jj);
%             valC(II) = lambda.*hCon(:,ii,jj);
%             n = n + 1;
%         end
%     end
%     
%     h = sparse([row(valX~=0);row(valY~=0);row(valR~=0);row(valC~=0)],...
%                [col(valX~=0);col(valY~=0);col(valR~=0);col(valC~=0)],...
%                [valX(valX~=0);valY(valY~=0);valR(valR~=0);valC(valC~=0)],...
%                4*NC,4*NC);
    
end
