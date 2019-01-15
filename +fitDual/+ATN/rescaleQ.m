function [ Q, iCells ] = rescaleQ( Q, iCells, Struct, extCell )
    % RESCALE Q 

    [ tri, bulk0, ext0, ~, ~, rV ] = fitDual.returnGraph( Struct, extCell );

    % Match iCells to our graph fitting ordering.
    idxs = arrayfun(@(x)find(iCells==x,1),[bulk0,ext0]);
    Q = Q(idxs,:);
    iCells = iCells(idxs);
    
    [ tri ] = fitDual.orderTri( Q, tri );
    
    % Subtract mean.
    R = mean(rV,1);
    rV = bsxfun(@minus,rV,R);
    
    % Find voronoi positions.
    [ r0 ] = fitDual.ATN.returnVertexPositions( Q, zeros(size(Q,1),1), tri );
    
    % Find scale between them.
    lambda = mean(dot(rV,r0,2))./mean(dot(r0,r0,2));
    
    Q = lambda*Q;
    Q = bsxfun(@plus,Q,R);
    
end

