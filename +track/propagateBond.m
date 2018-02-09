function [ bTrack ] = propagateBond( Struct, cPair, T )
    %PROPAGATE BOND 
        
    [ cTrack ] = track.propagateCells( cPair, T );
    bTrack = zeros(length(Struct(T(2)).Bdat),1);
    
    
    b0Cells = vertcat(Struct(T(1)).Bdat.cells); 
    bCells = vertcat(Struct(T(2)).Bdat.cells);
    bCells = cTrack(bCells);
    
    b0Cells = sort(b0Cells,2);
    bCells = sort(bCells,2);
    
    [~,ind0,ind] = intersect(b0Cells,bCells,'rows');
    bTrack(ind) = ind0;
    
end

