function [ hR ] = returnBondHessianTheta( q, theta, p, bCells )
    
    %% Calculculate elements of Rho and Radius Hessian
%     hRhoX = fitDual.AFN.deriv.rhoXHessTheta(p(bCells(:,1)),p(bCells(:,2)),q(bCells(:,1),1),q(bCells(:,2),1));
%     hRhoY = fitDual.AFN.deriv.rhoYHessTheta(p(bCells(:,1)),p(bCells(:,2)),q(bCells(:,1),2),q(bCells(:,2),2));
    hR = fitDual.AFN.deriv.radiusHessTheta(p(bCells(:,1)),p(bCells(:,2)),q(bCells(:,1),1), ...
                 q(bCells(:,2),1),q(bCells(:,1),2),q(bCells(:,2),2),theta(bCells(:,1)),theta(bCells(:,2)));
                         
%     hRhoX = reshape(hRhoX,size(hRhoX,1),2,2);
%     hRhoY = reshape(hRhoY,size(hRhoY,1),2,2);
    hR = reshape(hR,size(hR,1),2,2);

%     %% Add together in sparse matrix.
%     NB = size(bCells,1);
%     NC = size(q,1);
%     
%     row = zeros(3*NB,1);
%     col = zeros(3*NB,1);
% %     valX = zeros(3*NB,1);
% %     valY = zeros(3*NB,1);
%     valR = zeros(3*NB,1);
%     
%     n = 1;
%     for ii = 1:2
%         for jj = ii:2
%             
%             II = ((n-1)*NB) + (1:NB);
%             
%             if (ii <= 1)
%                 row(II) = bCells(:,1);
%             else
%                 row(II) = bCells(:,2);
%             end
%             
%             if (jj <= 1)
%                 col(II) = bCells(:,1);
%             else
%                 col(II) = bCells(:,2);
%             end
%             
% %             valX(II) = dNormX.*hRhoX(:,ii,jj);
% %             valY(II) = dNormY.*hRhoY(:,ii,jj);
%             valR(II) = -dAvg.*hR(:,ii,jj);
%             n = n + 1;
%         end
%     end
%     
%     h = sparse(row(valR~=0),col(valR~=0),valR(valR~=0),NC,NC);
    
end
