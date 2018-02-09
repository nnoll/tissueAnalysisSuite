function [ bInd, cInd ] = matchBondsAndCells( Struct, PN, extCell )
%MATCHBONDS Summary of this function goes here
%   Detailed explanation goes here

    [ dC, ~, ~, ~, bulkCells ] = fitDual.ATN.computeDiffOperators( Struct, extCell ); 
    
    rC = zeros(length(bulkCells),2);
    for ii = 1:length(bulkCells)
        rC(ii,:) = Struct.Cdat(bulkCells(ii)).centroid.coord;
    end
    xC = PN.returnCentroids();
    
    rB = .5*abs(dC)*rC;
    xB = .5*abs(PN.d0)*xC;
    
    Db = pdist2(xB,rB);
    bInd = track.munkres(Db);

    Dc = pdist2(xC,rC);
    cInd = track.munkres(Dc);
    
end

