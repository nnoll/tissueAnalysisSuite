function [ Z, numConComps, bulkCells ] = numberOfBulkCells( Struct, XRange, YRange )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    C = length(Struct.Cdat);
    bulkCells = 1:C;
    three_fold = [Struct.Cdat.all_threefold];
    if (isfield(Struct.Cdat,'hole'))
        holes = [Struct.Cdat.hole];
    else
        holes = zeros(size(bulkCells));
        holes(1) = 1;
    end

    if (isfield(Struct,'labelMat'))
        L = Struct.labelMat;
        YRange = round(YRange);
        XRange = round(XRange);
        containedCells = L(YRange(1):YRange(2),XRange(1):XRange(2));
        containedCells = unique(containedCells(containedCells>1));
    else
        rV = horzcat([Struct.Vdat.vertxcoord]',[Struct.Vdat.vertycoord]');
        containedVerts = ( rV(:,1) >= XRange(1) ) & ( rV(:,1) <= XRange(2) ) & ...
                         ( rV(:,2) >= YRange(1) ) & ( rV(:,2) <= YRange(2) );

        containedCells = [];
        for v = find(containedVerts)'
            containedCells = [containedCells,Struct.Vdat(v).ncells];
        end
        containedCells = unique(containedCells);
    end
    
    
    involvedCells = (three_fold==1) & (holes==0) & ismember(bulkCells,containedCells);  

    bndryCells = holes == 1;
    bndryCells = Struct.Cdat(bndryCells);
    bndryCells = unique([bndryCells.ncells]);
    
    bulkCells = bulkCells(involvedCells);
    bulkCells(bulkCells == 1) = [];
    bulkCells(ismember(bulkCells,bndryCells)) = [];
    bad_cells = [];
    for ii = 1:length(bulkCells)
        if (sum(ismember(Struct.Cdat(bulkCells(ii)).ncells,bulkCells))==0)
            bad_cells = horzcat(bad_cells,ii);
        end
    end
    bulkCells(bad_cells) = [];

    % Bulk verts are verts connected to bulk cells
    bulkVerts = [];
    for bc = bulkCells
        bulkVerts = horzcat(bulkVerts,Struct.Cdat(bc).nverts);
    end
    bulkVerts = unique(bulkVerts);

    % Find extCells
    involvedCells = [];
    for bv = bulkVerts
        involvedCells = horzcat(involvedCells,Struct.Vdat(bv).ncells);
    end
    involvedCells = unique(involvedCells);
    extCells = involvedCells(~ismember(involvedCells,bulkCells));
    involvedCells = [bulkCells,extCells];
    
    Z = length(involvedCells);
    if (Z >= 15)
        [ ~, dV, ~, ~, ~ ] = fitDual.subATN.computeSubDiffOperators( Struct, XRange, YRange );
        L = dV' * dV;
        numConComps = size(L,1) - rank(L);
    else
        numConComps = 0;
    end
   
end

