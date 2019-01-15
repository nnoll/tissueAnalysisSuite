function [ E, dE ] = energy( x, rV, tri, dBC, tBC, rBC )
    % ENERGY Attach springs between predicted vertex and observed vertex
    % positions.

    x = reshape(x,length(x)/4,4);
    q = x(:,1:2);
    theta = x(:,3);
    p = x(:,4);
    
    [ pTri ] = fitDual.AFN.returnPTri( q, p, tri );
    
    if (nargout == 2)
        [ r, gradR ] = fitDual.AFN.returnVertexPositions( q, theta, p, tri, pTri, rV );
    else
        [ r ] = fitDual.AFN.returnVertexPositions( q, theta, p, tri, pTri, rV );
    end
    
    E = mean( sum((r-rV).^2,2) )/2;

    if (nargin > 3)
       B = bsxfun(@times,dBC*p,rBC) - dBC*bsxfun(@times,p,q);
       E = E + .5*mean( (dot(B,tBC,2)).^2 ); 
    end
    
    if (nargout == 2)
        dE = zeros(size(gradR,1),4*size(q,1));

        NV = size(tri,1);
        NC = size(q,1);

        for ii = 1:3
            for jj = 1:4
                 dE( (1:NV) + NV*(NC*(jj-1)+pTri(:,ii)-1)' ) = gradR(:,ii,jj);
            end
        end

        dE = sum(dE,1)'/NV;
        
        if (nargin > 3)
            dE2 = zeros(size(B,1),size(q,1),4);
            dPos = dBC == 1;
            dNeg = dBC == -1;
            
            % Qx
                
            % Qy
            
            % Theta
            
            % p
            
            dE = dE + dE2;
        end
        
    end
    
end

