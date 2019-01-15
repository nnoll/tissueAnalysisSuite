function [ sigma, iCells ] = stressTensor( Struct, L, T )

    if (nargin == 2)
        [ T ] = MI.invertMech(Struct,1,1);
    end
    
    [ dC, dV, ~, iVerts, iCells ] = fitDual.ATN.computeDiffOperators( Struct, 1 );  
    rv = zeros(length(iVerts),2);
    for ii = 1:length(iVerts)
        rv(ii,1) = double(Struct.Vdat(iVerts(ii)).vertxcoord);
        rv(ii,2) = double(Struct.Vdat(iVerts(ii)).vertycoord);
    end
    
    rb = dV * rv;
    D = sqrt(sum(rb.^2,2));
    rb = bsxfun(@rdivide, rb, D);
    
    sigmaB = zeros(size(rb,1),3);
    sigmaB(:,1) = rb(:,1) .* T .* rb(:,1);
    sigmaB(:,2) = rb(:,1) .* T .* rb(:,2);
    sigmaB(:,3) = rb(:,2) .* T .* rb(:,2);

    sigma = sparse(abs(dC)') * sigmaB;
    S = regionprops(L,'Area');
    A = [S(iCells).Area];
    sigma = bsxfun(@rdivide,sigma,A');
    
end

