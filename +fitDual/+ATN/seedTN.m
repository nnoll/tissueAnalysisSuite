function [ QI ] = seedTN( Struct, extCell )
    % SEED TN 

    % Infer rough estimate of mechanics.
    T = MI.invertMech(Struct,extCell);
    
    % Build matrix for least squares fitting.
    Rot = [0,-1;1,0];
    [ dC, dV, ~, iVerts, iCells ] = fitDual.ATN.computeDiffOperators( Struct, extCell );
    
    rv = zeros(length(iVerts),2);
    for ii = 1:length(iVerts)
        rv(ii,1) = double(Struct.Vdat(iVerts(ii)).vertxcoord);
        rv(ii,2) = double(Struct.Vdat(iVerts(ii)).vertycoord);
    end

    %Obtain bond vectors.
    rb = dV * rv;
    D = sqrt(sum(rb.^2,2));
    nb = rb*Rot;
    nb = bsxfun(@rdivide,nb,D);
    rb = bsxfun(@rdivide,rb,D);
    
    rC =  zeros(length(iCells),2);
    for ii = 1:length(iCells)
        rC(ii,:) = Struct.Cdat(iCells(ii)).centroid.coord;
    end

    nC = dC * rC;
    S = sign(dot(nC,nb,2));
    dC = bsxfun(@times,S,dC);
        
    M = [bsxfun(@times,dC,nb(:,1)),bsxfun(@times,dC,nb(:,2))];
    M = [M;[bsxfun(@times,dC,rb(:,1)),bsxfun(@times,dC,rb(:,2))]];
    M = [M;[ones(1,size(dC,2)),zeros(1,size(dC,2))];[zeros(1,size(dC,2)),ones(1,size(dC,2))]];
    M = sparse(M);
    b = [T;zeros(size(T));0;0];
    QI = M \ b;
    QI = reshape(QI,length(QI)/2,2);

end

