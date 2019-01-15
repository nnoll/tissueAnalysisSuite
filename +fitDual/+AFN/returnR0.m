function [ r0, S0, etaStar, dCLin ] = returnR0( q, p, theta, pTri )
%RETURNR0 Summary of this function goes here
%   Detailed explanation goes here

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

    r0 =  (bsxfun(@times,(p1.*(sum(q1.^2,2)) + z1),t1) ... 
            + bsxfun(@times,(p2.*(sum(q2.^2,2)) + z2),t2) ...
            + bsxfun(@times,(p3.*(sum(q3.^2,2)) + z3),t3))*Rot';
    r0 = bsxfun(@rdivide,r0,4*S0);

    eta = bsxfun(@times,p1,t1) + bsxfun(@times,p2,t2) + bsxfun(@times,p3,t3);
    etaMag = sqrt(sum(eta.^2,2));
    etaStar = eta*Rot';
    
    if (nargout > 3) % Compute derivative components
        
        dCLin = zeros(size(r0,1),3,4);
        
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

            % Qx component
            dr0 = bsxfun(@times, p(a)./(4*S0), bsxfun(@times,dQ,Rot(:,1)') + 2*bsxfun(@times,q(a,1),tStar) - 2*bsxfun(@times,tStar(:,1),r0));
            dEta = bsxfun(@times,p(a).*(p(b) - p(c)),Rot(:,1)');
            dS = p(a).*tStar(:,1)/2;
            
            dCLin(:,ii,1) = 2 * (dot(etaStar,dEta,2) .* dot(r0,r0,2) + etaMag.^2 .*dot(r0,dr0,2) ...
                              - (2*S0-dot(etaStar,r0,2)).*(2*dS - dot(r0,dEta,2) - dot(etaStar,dr0,2)));

            % Qy component
            dr0 = bsxfun(@times, p(a)./(4*S0), bsxfun(@times,dQ,Rot(:,2)') + 2*bsxfun(@times,q(a,2),tStar) - 2*bsxfun(@times,tStar(:,2),r0));
            dEta = bsxfun(@times,p(a).*(p(b) - p(c)),Rot(:,2)');
            dS = p(a).*tStar(:,2)/2;
            
            dCLin(:,ii,2) = 2 * (dot(etaStar,dEta,2) .* dot(r0,r0,2) + etaMag.^2 .*dot(r0,dr0,2) ...
                              - (2*S0-dot(etaStar,r0,2)).*(2*dS - dot(r0,dEta,2) - dot(etaStar,dr0,2)));
                         
            % Theta component
            dr0 = bsxfun(@times, 1./(4*S0), tStar);
            
            dCLin(:,ii,3) = 2 * etaMag.^2.*dot(r0,dr0,2) + 2*(2*S0-dot(etaStar,r0,2)).*dot(etaStar,dr0,2);
            
            % P component
            dr0 = bsxfun(@times, 1./(4*S0), bsxfun(@times,Qsq,tStar) + bsxfun(@times,dQ,qStar) - 2*bsxfun(@times,dot(q(a,:),tStar,2),r0));
            dEta = tStar + bsxfun(@times,(p(b)-p(c)),qStar);
            dS = dot(q(a,:),tStar,2)/2;
            dCLin(:,ii,4) = 2 * (dot(etaStar,dEta,2) .* dot(r0,r0,2) + etaMag.^2 .*dot(r0,dr0,2) ...
                              - (2*S0-dot(etaStar,r0,2)).*(2*dS - dot(r0,dEta,2) - dot(etaStar,dr0,2)));
                         
        end
        
    end
    
end

