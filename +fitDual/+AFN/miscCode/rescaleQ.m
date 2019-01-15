function [ q ] = rescaleQ( q, p, Struct )
    % RESCALE Q 

    [ tri, ~, ~, ~, ~, rV ] = fitATN.returnGraph( Struct );
    pTri = fitAFN.returnPTri(q,p,tri);
    
    % Subtract mean.
    R = mean(rV,1);
    rV = bsxfun(@minus,rV,R);
    
    % Find voronoi positions.
    [ r0 ] = fitAFN.returnVertexPositions( q, zeros(size(p)), p, pTri, rV );
    
    % Find scale between them.
    lambda = mean(dot(rV,r0,2))./mean(dot(r0,r0,2));
    q = lambda*q;
    q = bsxfun(@plus,q,R);
    
end

