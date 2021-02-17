function [ E, dE ] = energy( x, rV, tri, A, b )
    % ENERGY Attach springs between predicted vertex and observed vertex
    % positions.

    x = reshape(x,length(x)/3,3);
    q = x(:,1:2);
    theta = x(:,3);
    
    [ tri ] = fitDual.orderTri(q,tri);
    [ r, gradR ] = fitDual.ATN.returnVertexPositions( q, theta, tri, rV );
   
    E = mean( sum((r-rV).^2,2) )/2;
    
    if (nargin == 5)
        E = E + .5 * mean((A*x(:)-b).^2);
    end
    
    if (nargout >= 2)
    % dE = zeros(size(gradR,1),3*size(q,1));

        NV = size(tri,1);
        NC = size(q,1);
        row = zeros(9*NV,1);
        col = zeros(9*NV,1);
        val = zeros(9*NV,1);

        n = 1;
        for ii = 1:3 % Iterate over cells
           for jj = 1:3 % Iterate over dual dof
                II = ((n-1)*NV) + (1:NV);
                row(II) = 1:NV;
                col(II) = (NC*(jj-1)) + tri(:,ii);
                val(II) = gradR(:,ii,jj)/NV;
                n = n + 1;
           end
        end
        
        dE = sparse(row,col,val,size(gradR,1),3*size(q,1));
%         for ii = 1:3
%             for jj = 1:3
%                  dE( (1:NV) + NV*(NC*(jj-1)+tri(:,ii)-1)' ) = gradR(:,ii,jj);
%             end
%         end

        dE = sum(dE,1)';
        
        if (nargin == 5)
            nCon = length(b);
            dE = dE + ((A*x(:) - b)'*A)'/nCon;
        end

    end

end

