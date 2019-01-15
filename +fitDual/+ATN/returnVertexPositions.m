function [ r, gradR ] = returnVertexPositions( q, theta, tri, r0 )
    % RETURN VERTEX POSITIONS 

    Rot = [0,-1;1,0];

    z1 = theta(tri(:,1));
    z2 = theta(tri(:,2));
    z3 = theta(tri(:,3));
    
    q1 = q(tri(:,1),:);
    q2 = q(tri(:,2),:);
    q3 = q(tri(:,3),:);
    
    t1 = q3 - q2;
    t2 = q1 - q3;
    t3 = q2 - q1;
    
    S0 = .5*abs(t1(:,1).*t2(:,2) - t1(:,2).*t2(:,1));
    
    r =  (bsxfun(@times,sum(q1.^2,2) + z1,t1) ... 
        + bsxfun(@times,sum(q2.^2,2) + z2,t2) ...
        + bsxfun(@times,sum(q3.^2,2) + z3,t3))*Rot';
    r = bsxfun(@rdivide,r,4*S0);
    
    if (nargout == 2)
        deltaR = r - r0;
        gradR = zeros(size(tri,1),3,3);   
        for ii = 1:3
            
            pII = mod(ii,3) + 1;
            mII = mod(ii-2,3) + 1;
            
            a = tri(:,ii);
            b = tri(:,pII);
            c = tri(:,mII);
            
            dQ = (sum(q(b,:).^2,2)+theta(b)) - (sum(q(c,:).^2,2)+theta(c));
            tStar = (q(c,:)-q(b,:))*Rot';
            
            dr0_x = bsxfun(@times, 1./(4*S0), bsxfun(@times,dQ,Rot(:,1)') + 2*bsxfun(@times,q(a,1),tStar) - 2*bsxfun(@times,tStar(:,1),r));
            dr0_y = bsxfun(@times, 1./(4*S0), bsxfun(@times,dQ,Rot(:,2)') + 2*bsxfun(@times,q(a,2),tStar) - 2*bsxfun(@times,tStar(:,2),r));
            dr0_t = bsxfun(@times, 1./(4*S0), tStar);
            
            gradR(:,ii,1) = dot(deltaR,dr0_x,2);
            gradR(:,ii,2) = dot(deltaR,dr0_y,2);
            gradR(:,ii,3) = dot(deltaR,dr0_t,2);

        end
        
    end
    
end

