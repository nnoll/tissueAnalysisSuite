function [ T, P ] = returnActualMech( Struct, extCell )
%RETURNACTUALMEH Summary of this function goes here
%   Detailed explanation goes here

    [ dC, dV, ~, iVerts, iCells ] = fitDual.ATN.computeDiffOperators( Struct, extCell );  
    T = zeros(size(dC,1),1);
    P = zeros(length(iCells),1);
    
    bVerts = horzcat(Struct.Bdat.verts)';
    for b = 1:size(dC,1)
        v1 = iVerts(dV(b,:)==1);
        v2 = iVerts(dV(b,:)==-1);
        B = find( ( (bVerts(:,1) == v1) & (bVerts(:,2) == v2) ) | ...
                  ( (bVerts(:,2) == v1) & (bVerts(:,1) == v2) ) );
        T(b) = Struct.Bdat(B).tension;
    end
    
    ii = 1;
    for c = iCells
       P(ii) = Struct.Cdat(c).p; 
       ii = ii + 1;
    end
    
end

