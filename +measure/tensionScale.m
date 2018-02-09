function [ scale ] = tensionScale( Struct, xG, yG )
    % TENSION SCALE 

    xB = [];
    T = [];
    for b = 1:length(Struct.Bdat)
        if (~isempty(Struct.Bdat(b).tension))
            T = [T;Struct.Bdat(b).tension];
            x = mean([Struct.Vdat(Struct.Bdat(b).verts).vertxcoord]);
            y = mean([Struct.Vdat(Struct.Bdat(b).verts).vertycoord]);
            
            xB = [xB;[x,y]];
        end
    end

    scale = zeros(size(yG,2),size(xG,2));
    for ii = 1:size(yG,2)
        for jj = 1:size(xG,2)
            bonds = (xB(:,1) > xG(1,jj)) & (xB(:,1) < xG(2,jj)) & (xB(:,2) > yG(1,ii)) & (xB(:,2) < yG(2,ii));
            scale(ii,jj) = nanmedian(T(bonds));
        end
    end
end

