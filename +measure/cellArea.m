function [ A ] = cellArea( Struct, L )
    % CELL AREA
    
    for t = 1:length(Struct)
       S = regionprops(L,'Area');
       A{t} = [S.Area];
    end

end

