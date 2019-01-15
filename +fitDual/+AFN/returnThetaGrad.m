function [ dR, dR0, r0, etaStar ] = returnThetaGrad( r, q, p, pTri, indB, r0, lambda )
% RETURNTHETAGRAD 

    if (nargin > 4)
        mode = 0;
    else
        mode = 1;
    end
    
    % Compute eta;
    Rot = [0,-1;1,0];
            
    p1 = p(pTri(:,1));
    p2 = p(pTri(:,2));
    p3 = p(pTri(:,3));

    q1 = q(pTri(:,1),:);
    q2 = q(pTri(:,2),:);
    q3 = q(pTri(:,3),:);

    % Store triangular edges;
    t1 = bsxfun(@times,p3,q3) - bsxfun(@times,p2,q2);
    t2 = bsxfun(@times,p1,q1) - bsxfun(@times,p3,q3);
    t3 = bsxfun(@times,p2,q2) - bsxfun(@times,p1,q1);

    S0 = .5*abs(t1(:,1).*t2(:,2) - t1(:,2).*t2(:,1));

    eta = bsxfun(@times,p1,t1) + bsxfun(@times,p2,t2) + bsxfun(@times,p3,t3);
    etaMag = sqrt(sum(eta.^2,2));
    etaStar = eta*Rot';
    
    if (nargin > 4)
        Denom = 2*S0 - (lambda .* etaMag.^2) - dot(r0,etaStar,2);
    end
    
    dR = zeros(size(r,1),2,3);
    dR0 = zeros(size(r,1),2,3);

    for ii = 1:3
        
       pII = mod(ii,3) + 1;
       mII = mod(ii-2,3) + 1;

       a = pTri(:,ii);
       b = pTri(:,pII);
       c = pTri(:,mII); 
       tStar = (bsxfun(@times,p(c),q(c,:))-bsxfun(@times,p(b),q(b,:)))*Rot';

       % r0 component
       dR0(:,:,ii) = bsxfun(@times, 1./(4*S0), tStar);

       if (mode == 0)
            dr0_t = dR0(:,:,ii);

            gradLambda = dot(r,dr0_t,2) ./ Denom;
            altGrad = -dot(etaStar,dr0_t,2) ./ etaMag.^2;
            gradLambda(indB) = altGrad(indB);
            dR(:,:,ii) = dR0(:,:,ii) + bsxfun(@times,etaStar,gradLambda);
            
       end
       
    end
    
end

