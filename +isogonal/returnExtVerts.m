function [ extVerts, rVext ] = returnExtVerts( bulkVerts, Struct )
    % RETURN EXT VERTS 

    iVerts = zeros(3*length(bulkVerts),1);
    for ii = 1:length(bulkVerts)
        iVerts( (3*(ii-1) + 1):(3*ii)  ) = Struct.Vdat(bulkVerts(ii)).nverts;
    end
    
    iVerts = unique(iVerts);
    extVerts = iVerts(~ismember(iVerts,bulkVerts))';
    
    rVext = zeros(length(extVerts),2);
    for ii = 1:length(extVerts)
        rVext(ii,:) = [double(Struct.Vdat(extVerts(ii)).vertxcoord),double(Struct.Vdat(extVerts(ii)).vertycoord)];
    end
    
end

