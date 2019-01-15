function [ Z ] = numberOfBulkCells( Struct, XRange, YRange )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

    C = length(Struct.Cdat);
    bulkCells = 1:C;
    three_fold = [Struct.Cdat.all_threefold];
    holes = [Struct.Cdat.hole];

%     Rc = [Struct.Cdat.centroid];
%     Rc = vertcat(Rc.coord);
%     involvedCells = (three_fold==1) & (Rc(:,1)' >= XRange(1)) & (Rc(:,1)' <= XRange(2)) & ...
%                     (Rc(:,2)' >= YRange(1)) & (Rc(:,2)' <= YRange(2));
    YRange = round(YRange);
    XRange = round(XRange);
    L = Struct.labelMat;
    containedCells = L(YRange(1):YRange(2),XRange(1):XRange(2));
    containedCells = unique(containedCells(containedCells>1));
    
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

    Z = length(bulkCells);
    
end

