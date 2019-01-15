function [ Struct, i1 ] = storeMech( Struct, extCell, mode )
%STOREMECH Summary of this function goes here
%   Detailed explanation goes here

    for t = 1:length(Struct)
        t
        [ T, P ] = MI.invertMech( Struct(t), extCell, mode );
        [ dC, ~, ~, ~, iCells ] = fitDual.ATN.computeDiffOperators( Struct(t), extCell );  
        i1{t} = zeros(length(Struct(t).Bdat),1);
        for b = 1:length(Struct(t).Bdat)
            bCells = Struct(t).Bdat(b).cells;
            if (sum(ismember(bCells,iCells)) == 2) % Both cells involved in fit.
                cellInd1 = (ismember(iCells,bCells(1)));
                cellInd2 = (ismember(iCells,bCells(2)));
                bInd = ( dC(:,cellInd1)~=0 &  dC(:,cellInd2)~=0 );
                Struct(t).Bdat(b).tension = T(bInd);
                i1{t}(b) = 1;
            end
        end
        
        for c = 1:length(iCells)
            Struct(t).Cdat(iCells(c)).pressure = P(c);
        end
        i1{t} = find(i1{t});
    end

end

