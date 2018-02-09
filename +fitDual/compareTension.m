function [ T, Ta ] = compareTension( Struct, PN )
%COMPARETENSION Summary of this function goes here
%   Detailed explanation goes here

    [ d0, ~, ~, ~, i0 ] = fitDual.ATN.computeDiffOperators( Struct, 1 );
    [ bCells ] = generate.bCells( d0, i0 );
    [ i1 ] = generate.bondMap( Struct );

    T =  PN.returnTension;
    [ bGCells ] = generate.bCells( PN.d0, PN.cellLabels );
    bInd = zeros(size(bGCells,1),1);
    for b = 1:size(bGCells,1)
        bInd(b) = find(ismember(bCells,bGCells(b,:),'rows'));
    end
    Ta = [Struct.Bdat(i1{1}(bInd)).actual_tension]';
    
end

