function [ r, indB, r0, lambda, gradR ] = returnVertexPositions( q, theta, p, tri, pTri, rV )
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

    r0 =  (bsxfun(@times,p1.*sum(q1.^2,2) + z1,t1) ... 
         + bsxfun(@times,p2.*sum(q2.^2,2) + z2,t2) ...
         + bsxfun(@times,p3.*sum(q3.^2,2) + z3,t3))*Rot';
    r0 = bsxfun(@rdivide,r0,4*S0);
    
    % Compute the deviation off the zeroth order term.
    eta = bsxfun(@times,p1,t1) + bsxfun(@times,p2,t2) + bsxfun(@times,p3,t3);
    etaMag = sqrt(sum(eta.^2,2));

    etaStar = eta*Rot';
    etaStar = bsxfun(@rdivide,etaStar,etaMag);
    
    A = ((2*S0./etaMag) - dot(etaStar,r0,2));
    B = A.^2 - dot(r0,r0,2); %.*dot(etaStar,etaStar,2);
    indB = B<0;
    B(indB) = 0;

    lambdaMinus = A - sqrt( B );
    lambdaPlus = A + sqrt( B );
    
    lambdaMinus = real(lambdaMinus);
    lambdaPlus = real(lambdaPlus);
    
%     ind3 = isnan(etaStar);
%     etaStar(ind3) = 0;
    
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
    lambda = bsxfun(@rdivide,lambda,etaMag);
    if (nargout == 5)
        deltaR = r - rV;
        gradR = zeros(size(pTri,1),3,4);
        etaStar = bsxfun(@times,etaStar,etaMag);
        for ii = 1:3
            
            pII = mod(ii,3) + 1;
            mII = mod(ii-2,3) + 1;
            
            a = pTri(:,ii);
            b = pTri(:,pII);
            c = pTri(:,mII);
            
            dQ = p(b).*(sum(q(b,:).^2,2))+theta(b) - p(c).*(sum(q(c,:).^2,2))+theta(c);
            Qsq = (sum(q(a,:).^2,2));
            qStar = q(a,:)*Rot';
            tStar = (bsxfun(@times,p(c),q(c,:))-bsxfun(@times,p(b),q(b,:)))*Rot';
            
            %% Contribution from the zeroth order term.
            
            % dQx
            dr0_x = bsxfun(@times, p(a)./(4*S0), bsxfun(@times,dQ,Rot(:,1)') + 2*bsxfun(@times,q(a,1),tStar) - 2*bsxfun(@times,tStar(:,1),r0));
            dEta_x = bsxfun(@times,p(a).*(p(b) - p(c)),Rot(:,1)');
            dS_x = p(a).*tStar(:,1)/2;

            gradR(:,ii,1) = dot(deltaR,dr0_x,2);
                  
            % dQY
            dr0_y = bsxfun(@times, p(a)./(4*S0), bsxfun(@times,dQ,Rot(:,2)') + 2*bsxfun(@times,q(a,2),tStar) - 2*bsxfun(@times,tStar(:,2),r0));
            dEta_y = bsxfun(@times,p(a).*(p(b) - p(c)),Rot(:,2)');
            dS_y = p(a).*tStar(:,2)/2;

            gradR(:,ii,2) = dot(deltaR,dr0_y,2);
                   
            % dTheta
            dr0_t = bsxfun(@times, 1./(4*S0), tStar);
            dEta_t = zeros(size(dr0_t));
            dS_t = zeros(size(dS_y));

            gradR(:,ii,3) = dot(deltaR,dr0_t,2);
            
            % dP
            dr0_p = bsxfun(@times, 1./(4*S0), bsxfun(@times,Qsq,tStar) + bsxfun(@times,dQ,qStar) - 2*bsxfun(@times,dot(q(a,:),tStar,2),r0));
            dEta_p = tStar + bsxfun(@times,(p(b)-p(c)),qStar);
            dS_p = dot(q(a,:),tStar,2)/2;
            
            gradR(:,ii,4) = dot(deltaR,dr0_p,2);

            % Contribution from the first term of the delta term.
            
            gradR(:,ii,1) = gradR(:,ii,1) + lambda .* dot(dEta_x,deltaR,2);
            gradR(:,ii,2) = gradR(:,ii,2) + lambda .* dot(dEta_y,deltaR,2);
            gradR(:,ii,3) = gradR(:,ii,3) + lambda .* dot(dEta_t,deltaR,2);
            gradR(:,ii,4) = gradR(:,ii,4) + lambda .* dot(dEta_p,deltaR,2);
            
            % Contribution from the second term of the delta term.
            Denom = 2*S0 - (lambda .* etaMag.^2) - dot(r0,etaStar,2);
            Indx = B == 0;
            edgeCase = any(Indx);
            
            % dQx
            term1 = (dot(r,dr0_x + bsxfun(@times,lambda,dEta_x),2) - 2*lambda.*dS_x)./Denom;
            if (edgeCase)
                tmp = (2*dS_x - dot( bsxfun(@times,2*lambda,etaStar) + r0, dEta_x, 2) - dot(etaStar,dr0_x,2)) ./ etaMag.^2;
                term1(Indx) = tmp(Indx);
            end
            gradR(:,ii,1) = gradR(:,ii,1) + dot(deltaR,etaStar,2).*term1;
                  
            % dQY
            term2 = (dot(r,dr0_y + bsxfun(@times,lambda,dEta_y),2) - 2*lambda.*dS_y)./Denom;
            if (edgeCase)
                tmp = (2*dS_y - dot( bsxfun(@times,2*lambda,etaStar) + r0, dEta_y, 2) - dot(etaStar,dr0_y,2)) ./ etaMag.^2;
                term2(Indx) = tmp(Indx);
            end
            gradR(:,ii,2) = gradR(:,ii,2) + dot(deltaR,etaStar,2).*term2;
                   
            % dTheta
            term3 = (dot(r,dr0_t + bsxfun(@times,lambda,dEta_t),2) - 2*lambda.*dS_t)./Denom;
            if (edgeCase)
                tmp = (2*dS_t - dot( bsxfun(@times,2*lambda,etaStar) + r0, dEta_t, 2) - dot(etaStar,dr0_t,2)) ./ etaMag.^2;
                term3(Indx) = tmp(Indx);
            end
            gradR(:,ii,3) = gradR(:,ii,3) + dot(etaStar,deltaR,2).*term3;
            
            % dP
            term4 = (dot(r,dr0_p + bsxfun(@times,lambda,dEta_p),2) - 2*lambda.*dS_p)./Denom;
            if (edgeCase)
                tmp = (2*dS_p - dot( bsxfun(@times,2*lambda,etaStar) + r0, dEta_p, 2) - dot(etaStar,dr0_p,2)) ./ etaMag.^2;
                term4(Indx) = tmp(Indx);
            end
            gradR(:,ii,4) = gradR(:,ii,4) + dot(deltaR,etaStar,2).*term4;
            
        end
    end     
    
end

