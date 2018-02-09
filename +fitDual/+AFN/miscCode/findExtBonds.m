function [ extBonds, rBext, rCext ] = findExtBonds( Struct, extCells )
    %FIND EXT BONDS 
    
    extBonds = [];
    rBext = [];
    rCext = [];
    for c = extCells
        nCells = Struct.Cdat(c).ncells;
        nExtCells = nCells(ismember(nCells,extCells));
        nExtCells = nExtCells(nExtCells > c);
        for nc = nExtCells
            bondVerts = Struct.Cdat(c).nverts(ismember(Struct.Cdat(c).nverts,Struct.Cdat(nc).nverts));
            
            if (length(bondVerts) == 2)
                extBonds = [extBonds;[find(extCells==c),find(extCells==nc)]];
                rBext = [rBext;[Struct.Vdat(bondVerts(1)).vertxcoord-Struct.Vdat(bondVerts(2)).vertxcoord, ...
                                Struct.Vdat(bondVerts(1)).vertycoord-Struct.Vdat(bondVerts(2)).vertycoord]];
                rCext = [rCext;Struct.Cdat(c).centroid.coord-Struct.Cdat(nc).centroid.coord];  
            end
        end
    end
    
    rBext = bsxfun(@rdivide,rBext,sqrt(sum(rBext.^2,2)).*sqrt(sum(rCext.^2,2)));
    
end

