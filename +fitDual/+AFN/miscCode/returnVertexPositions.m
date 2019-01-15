function [ r, gradR ] = returnVertexPositions( q, theta, p, tri, pTri, rV )
    % RETURN VERTEX POSITIONS 

    Rot = [0,-1;1,0];
    
    % Store pressure, theta, and q positions.
    p1 = p(pTri(:,1));
    p2 = p(pTri(:,2));
    p3 = p(pTri(:,3));
    
    z1 = theta(pTri(:,1));
    z2 = theta(pTri(:,2));
    z3 = theta(pTri(:,3));
    
    q1 = q(pTri(:,1),:);
    q2 = q(pTri(:,2),:);
    q3 = q(pTri(:,3),:);

    % Store triangular edges and its area.
    t1 = bsxfun(@times,p3,q3) - bsxfun(@times,p2,q2);
    t2 = bsxfun(@times,p1,q1) - bsxfun(@times,p3,q3);
    t3 = bsxfun(@times,p2,q2) - bsxfun(@times,p1,q1);
    
    S0 = .5*abs(t1(:,1).*t2(:,2) - t1(:,2).*t2(:,1));

    r0 =  (bsxfun(@times,(p1.*(sum(q1.^2,2) + z1)),t1) ... 
         + bsxfun(@times,(p2.*(sum(q2.^2,2) + z2)),t2) ...
         + bsxfun(@times,(p3.*(sum(q3.^2,2) + z3)),t3))*Rot';
    r0 = bsxfun(@rdivide,r0,4*S0);
    
    % Compute the deviation off the zeroth order term.
    eta = bsxfun(@times,p1,t1) + bsxfun(@times,p2,t2) + bsxfun(@times,p3,t3);
    etaMag = sqrt(sum(eta.^2,2));

    etaStar = eta*Rot';
    etaStar = bsxfun(@rdivide,etaStar,etaMag);
    
    A = ((2*S0./etaMag) - dot(etaStar,r0,2));
    B = A.^2 - dot(r0,r0,2); %.*dot(etaStar,etaStar,2);
    B(B<0) = 0;

    lambdaMinus = A - sqrt( B );
    lambdaPlus = A + sqrt( B );
    
%     lambdaMinus = lambdaMinus./etaMag.^2;
%     lambdaPlus = lambdaPlus./etaMag.^2;

%     ind1 = isnan(lambdaMinus);
%     lambdaMinus(ind1) = 0;
    lambdaMinus = real(lambdaMinus);
%     
%     ind2 = isnan(lambdaPlus);
%     lambdaPlus(ind2) = 0;
    lambdaPlus = real(lambdaPlus);
    
    ind3 = isnan(etaStar);
    etaStar(ind3) = 0;
    
    deltaR_M = bsxfun(@times,lambdaMinus,etaStar);
    deltaR_P = bsxfun(@times,lambdaPlus,etaStar);
        
    r1 = r0 + deltaR_M;
    r2 = r0 + deltaR_P;

    if (nargin == 6)
        D1 = sum((r1-rV).^2,2);
        D2 = sum((r2-rV).^2,2);
        [~,ind] = min([D1,D2],[],2);
        
        r = r1;
        r(ind==2,:) = r2(ind==2,:);  
        lambda = lambdaMinus;
        lambda(ind==2) = lambdaPlus(ind==2);
    else
        ind = sum(tri == pTri,2);
        ind = (ind == 1);
        deltaR_M(ind) = deltaR_P(ind);
    
        deltaR = bsxfun(@times,deltaR_M,etaStar);
        r = r0 + deltaR;
    end

    % Update vertex position as well as lambda for gradient calculation.

    if (nargout == 2)
        
        deltaR = r - rV;
        gradR = zeros(size(pTri,1),3,4);
        etaStar = bsxfun(@times,etaStar,etaMag);
        lambda = bsxfun(@rdivide,lambda,etaMag);
        
        for ii = 1:3
            
            pII = mod(ii,3) + 1;
            mII = mod(ii-2,3) + 1;
            
            a = pTri(:,ii);
            b = pTri(:,pII);
            c = pTri(:,mII);
            
            dQ = p(b).*(sum(q(b,:).^2,2)+theta(b)) - p(c).*(sum(q(c,:).^2,2)+theta(c));
            Qsq = (sum(q(a,:).^2,2)+theta(a));
            qStar = q(a,:)*Rot';
            tStar = (bsxfun(@times,p(c),q(c,:))-bsxfun(@times,p(b),q(b,:)))*Rot';
            
            %% Contribution from the zeroth order term.
            
            % dQx
            dr0_x = bsxfun(@times, p(a)./(4*S0), bsxfun(@times,dQ,Rot(:,1)') + 2*bsxfun(@times,q(a,1),tStar) - 2*bsxfun(@times,tStar(:,1),r0));
            dEta_x = bsxfun(@times,p(a).*(p(b) - p(c)),Rot(:,1)');
            dS_x = p(a).*tStar(:,1)/2;
            gradR(:,ii,1) = bsxfun(@times,dQ,deltaR(:,2));
            gradR(:,ii,1) = gradR(:,ii,1) + 2.*dot(tStar,deltaR,2).*q(a,1);
            gradR(:,ii,1) = gradR(:,ii,1) - 2.*dot(deltaR,r0,2).*tStar(:,1);
            gradR(:,ii,1) = p(a).*gradR(:,ii,1)./(4*S0);
                  
            % dQY
            gradR(:,ii,2) = bsxfun(@times,dQ,-deltaR(:,1));
            gradR(:,ii,2) = gradR(:,ii,2) + 2*dot(tStar,deltaR,2).*q(a,2);
            gradR(:,ii,2) = gradR(:,ii,2) - 2*dot(deltaR,r0,2).*tStar(:,2);
            gradR(:,ii,2) = p(a).*gradR(:,ii,2)./(4*S0);
                   
            % dTheta
            gradR(:,ii,3) = p(a) .* dot(deltaR,tStar,2)./(4*S0);
            
            % dP
            gradR(:,ii,4) = Qsq.*dot(tStar,deltaR,2);
            gradR(:,ii,4) = gradR(:,ii,4) + dQ.*dot(qStar,deltaR,2);
            gradR(:,ii,4) = gradR(:,ii,4) - 2*dot(r0,deltaR,2).*dot(q(a,:),tStar,2);
            gradR(:,ii,4) = gradR(:,ii,4)./(4*S0);
            
            %% Contribution from the first term of the delta term.
            
            % dQx
            gradR(:,ii,1) = gradR(:,ii,1) + lambda .* deltaR(:,2) .* p(a) .* (p(b)-p(c));
                  
            % dQY
            gradR(:,ii,2) = gradR(:,ii,2) - lambda .* deltaR(:,1) .* p(a) .* (p(b)-p(c));
                   
            % dTheta
            gradR(:,ii,3) = gradR(:,ii,3) + 0;
            
            % dP
            gradR(:,ii,4) = gradR(:,ii,4) + lambda .* (dot(deltaR,tStar,2) + dot(deltaR,bsxfun(@times,p(b)-p(c),qStar),2));
            
            %% Contribution from the second term of the delta term.
            Denom = 2*S0 - (lambda .* etaMag.^2) - dot(r0,etaStar,2);
            Indx = B == 0;
            edgeCase = any(Indx);
            
            % dQx
            term1 = lambda .* r(:,2) .* p(a) .* (p(b)-p(c));
            term1 = term1 + p(a) .* ( r(:,2).*dQ ) ./ (4*S0);
            term1 = term1 + p(a) .* ( q(a,1) .* dot(tStar,r,2) ) ./ (2*S0);
            term1 = term1 - p(a) .* ( tStar(:,1) .* dot(r0,r,2) ) ./ (2*S0);
            term1 = term1 - lambda .* p(a) .* tStar(:,1);
            term1 = term1 .* dot(etaStar,deltaR,2);
            term1 = term1./Denom;
            if (edgeCase)
                tmp = p(a) .* tStar(:,1);
                tmp = tmp - (2.*lambda.*etaStar(:,2) + r0(:,2)).*p(a).*(p(b)-p(c));
                tmp = tmp - p(a).*(dQ.*etaStar(:,2)+2*(q(a,1).*dot(etaStar,tStar,2)-tStar(:,1).*dot(r0,etaStar,2)))./(4*S0);
                tmp = (tmp .* dot(etaStar,deltaR,2)) ./ (etaMag.^2);
                term1(Indx) = tmp(Indx);
            end
            gradR(:,ii,1) = gradR(:,ii,1) + term1;
                  
            % dQY
            term2 = -lambda .* r(:,1) .* p(a) .* (p(b)-p(c)); 
            term2 = term2 - p(a) .* ( r(:,1).*dQ ) ./ (4*S0);
            term2 = term2 + p(a) .* ( q(a,2) .* dot(tStar,r,2) ) ./ (2*S0);
            term2 = term2 - p(a) .* ( tStar(:,2) .* dot(r0,r,2) ) ./ (2*S0);
            term2 = term2 - lambda .* p(a) .* tStar(:,2);
            term2 = term2 .* dot(etaStar,deltaR,2);
            term2 = term2./Denom;
            if (edgeCase)
                tmp = p(a) .* tStar(:,2);
                tmp = tmp + (2.*lambda.*etaStar(:,1) + r0(:,1)).*p(a).*(p(b)-p(c));
                tmp = tmp - p(a).*(-dQ.*etaStar(:,1)+2*(q(a,2).*dot(etaStar,tStar,2)-tStar(:,2).*dot(r0,etaStar,2)))./(4*S0);
                tmp = (tmp.* dot(etaStar,deltaR,2)) ./ (etaMag.^2);
                term2(Indx) = tmp(Indx);
            end
            gradR(:,ii,2) = gradR(:,ii,2) + term2;
                   
            % dTheta
            term3 = p(a) .* dot(tStar,r,2) ./ (4*S0);
            term3 = term3 .* dot(etaStar,deltaR,2);
            term3 = term3./Denom;
            if (edgeCase)
                tmp = -(p(a) .* dot(etaStar,tStar,2) .* dot(etaStar,deltaR,2)) ./ (4*S0.*etaMag.^2);
                term3(Indx) = tmp(Indx);
            end
            gradR(:,ii,3) = gradR(:,ii,3) + term3;
            
            % dP
            term4 = lambda .* (dot(r,tStar,2) + dot(r,bsxfun(@times,p(b)-p(c),qStar),2));
            
            r0Vec = bsxfun(@times,Qsq,tStar) + bsxfun(@times,dQ,qStar) - 2*bsxfun(@times,dot(q(a,:),tStar,2),r0);
            r0Vec = bsxfun(@rdivide,r0Vec,4*S0);
            
            term4 = term4 + dot(r,r0Vec,2);
            term4 = term4 - lambda .* dot(q(a,:),tStar,2);
            term4 = term4 .* dot(etaStar,deltaR,2);
            term4 = term4./Denom;
            if (edgeCase)
                tmp = dot(q(a,:),tStar,2);
                tmp = tmp - dot(bsxfun(@times,2*lambda,etaStar)+r0,tStar+bsxfun(@times,p(b)-p(c),qStar),2);
                tmp = tmp - (Qsq .* dot(tStar,etaStar,2) + dQ.*dot(etaStar,qStar,2) - 2*dot(r0,etaStar,2).*dot(q(a,:),tStar,2))./(4*S0);
                tmp = (tmp.* dot(etaStar,deltaR,2)) ./ (etaMag.^2);
                term4(Indx) = tmp(Indx);
            end
            gradR(:,ii,4) = gradR(:,ii,4) + term4;
            
        end
    end     
    
end

