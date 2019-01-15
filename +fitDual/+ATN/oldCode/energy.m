function [ E ] = energy( x, rV, tri, Aeq )
    % ENERGY Attach springs between predicted vertex and observed vertex
    % positions.

    x = reshape(x,length(x)/3,3);
    q = x(:,1:2);
    theta = x(:,3);
    
    [ tri ] = fitDual.orderTri(q,tri);
    [ r ] = fitDual.ATN.returnVertexPositions( q, theta, tri );
   
    if (nargin == 3)
        E = mean( sum((r-rV).^2,2) )/2;
    else
        E = .5*( mean( sum((r-rV).^2,2) ) + mean((Aeq*x(:)).^2) );
    end
    
end

