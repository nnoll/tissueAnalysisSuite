function [ Struct, i1 ] = uploadMechanics( this, Struct )
%UPLOADTENSIONS Summary of this function goes here
%   Detailed explanation goes here

    [ T ] = this.returnTension();
    [ p ] = this.p;
    
    i1 = zeros(length(Struct.Bdat),1);
    for b = 1:length(Struct.Bdat)
        bCells = Struct.Bdat(b).cells;
        if (sum(ismember(bCells,this.cellLabels)) == 2) % Both cells involved in fit.
            cellInd1 = find(ismember(this.cellLabels,bCells(1)));
            cellInd2 = find(ismember(this.cellLabels,bCells(2)));
            bInd = find( this.d0(:,cellInd1)~=0 &  this.d0(:,cellInd2)~=0 );
            if (isfield(Struct.Bdat(b),'chem'))
                if (~isempty(bInd) && ~isempty(Struct.Bdat(b).chem))
                    Struct.Bdat(b).tension = T(bInd);
                    i1(b) = 1;
                end
            else
                if (~isempty(bInd))
                    Struct.Bdat(b).tension = T(bInd);
                    i1(b) = 1;
                end
            end
        end
    end
    i1 = find(i1);
    
    for c = 1:length(p)
        cellLabel = this.cellLabels(c);
        Struct.Cdat(cellLabel).pressure = p(c);
    end
    
end

