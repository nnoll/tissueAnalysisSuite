function [ rB, rTh, grad_rB, grad_rTh ] = returnPartitionedVertexPositions( q, theta, tri )
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
    
    rB =  (bsxfun(@times,sum(q1.^2,2),t1) ... 
         + bsxfun(@times,sum(q2.^2,2),t2) ...
         + bsxfun(@times,sum(q3.^2,2),t3))*Rot';
    rB = bsxfun(@rdivide,rB,4*S0);
    
    rTh = (bsxfun(@times,z1,t1) ... 
         + bsxfun(@times,z2,t2) ...
         + bsxfun(@times,z3,t3))*Rot';
    rTh = bsxfun(@rdivide,rTh,4*S0);
    
    if (nargout > 2)
        
        grad_rB = zeros(size(tri,1),2,3,3);   
        grad_rTh = zeros(size(tri,1),2,3,3);   

        for ii = 1:3
            
            pII = mod(ii,3) + 1;
            mII = mod(ii-2,3) + 1;
            
            a = tri(:,ii);
            b = tri(:,pII);
            c = tri(:,mII);
            
            dQ = sum(q(b,:).^2,2) - sum(q(c,:).^2,2);
            dT = theta(b) - theta(c);
            tStar = (q(c,:)-q(b,:))*Rot';
            
            drB_x = bsxfun(@times, 1./(4*S0), bsxfun(@times,dQ,Rot(:,1)') + 2*bsxfun(@times,q(a,1),tStar) - 2*bsxfun(@times,tStar(:,1),rB));
            drB_y = bsxfun(@times, 1./(4*S0), bsxfun(@times,dQ,Rot(:,2)') + 2*bsxfun(@times,q(a,2),tStar) - 2*bsxfun(@times,tStar(:,2),rB));
            drB_t = zeros(size(drB_x));
            
            drTh_x = bsxfun(@times, 1./(4*S0), bsxfun(@times,dT,Rot(:,1)') - 2*bsxfun(@times,tStar(:,1),rTh));
            drTh_y = bsxfun(@times, 1./(4*S0), bsxfun(@times,dT,Rot(:,2)') - 2*bsxfun(@times,tStar(:,2),rTh));
            drTh_t = bsxfun(@times, 1./(4*S0), tStar);
            
            grad_rB(:,:,ii,1) = drB_x;
            grad_rB(:,:,ii,2) = drB_y;
            grad_rB(:,:,ii,3) = drB_t;
            
            grad_rTh(:,:,ii,1) = drTh_x;
            grad_rTh(:,:,ii,2) = drTh_y;
            grad_rTH(:,:,ii,3) = drTh_t;

        end
        
    end
    
end

