function [ clin, ceq, dlin, deq ] = thetaCon( x, q, p, pTri )
    
end

% [ r0, S0, etaStar ] = fitDual.AFN.returnR0( q, p, x, pTri );
%     
%     etaMag = sqrt(sum(etaStar.^2,2));
%     r0Mag = sqrt(dot(r0,r0,2));
%     clin = etaMag.*r0Mag + dot(etaStar,r0,2) - 2*S0;    
%     ceq = [];
%     
%     if (nargout == 4)
%         [ dR ] = fitDual.AFN.returnThetaGrad( r0, q, p, pTri );
%         dClin = zeros(size(r0,1),size(q,1));
%         NV = size(pTri,1);
%         for ii = 1:3
%             IPr = bsxfun(@times,r0(:,1),dR(:,1,ii)) + bsxfun(@times,r0(:,2),dR(:,2,ii));
%             IPe = bsxfun(@times,etaStar(:,1),dR(:,1,ii)) + bsxfun(@times,etaStar(:,2),dR(:,2,ii));
% 
%             dClin( (1:NV) + NV*(pTri(:,ii)-1)' ) = (etaMag .* IPr ./ r0Mag) + IPe; 
%         end
%         
%         dlin = sparse(dClin');
%         deq = [];
%         
%     end