function [ i1 ] = bondMap( Struct )
%BONDMAP

    i1 = cell(length(Struct),1);
    for t = 1:length(Struct)
        [ d0, ~, ~, ~, i0 ] = fitDual.ATN.computeDiffOperators( Struct(t), 1 );
        i1{t} = zeros(size(d0,1),1);

        for b = 1:length(Struct.Bdat)
            bCells = Struct.Bdat(b).cells;
            if (sum(ismember(bCells,i0)) == 2) % Both cells involved in diff operator.

                cellInd1 = ismember(i0,bCells(1));
                cellInd2 = ismember(i0,bCells(2));
                bInd = ( d0(:,cellInd1)~=0 &  d0(:,cellInd2)~=0 );

                i1{t}(bInd) = b;
            end

        end
    end


end

