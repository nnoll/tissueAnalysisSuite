function [ goodCells ] = goodCells( Struct, ERes, xG, yG )
    % GOOD CELLS 

    goodCells = 1:length(Struct.Cdat);
    badCells = [];
    [badRow,badCol] = find(ERes > 10);
%     imagesc(ERes)
%     hold on
%     scatter(badCol,badRow,'g','filled')
%     pause
    for n = 1:length(badRow)
        [ ~, ~, i0 ] = measure.numberOfBulkCells( Struct, xG(:,badCol(n)), yG(:,badRow(n)) );
        badCells = [badCells,i0];
    end
    
    badCells = [badCells,find([Struct.Cdat.hole] == 1)];
    badCells = unique(badCells);
    goodCells(ismember(goodCells,badCells)) = [];
    
end

