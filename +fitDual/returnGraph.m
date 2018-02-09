function [ Tri, bulkCells, extCells, bulkVerts, rC, rV ] = returnGraph( Struct, mode )
    
    if (nargin == 1)
        mode = 0;
    end
    
    if (mode == 0) % Don't exclude cell 1 
        
        % Find all bulk verts (three-fold).
        goodVerts = zeros(1,length(Struct.Vdat));
        for v = 1:length(goodVerts)
            goodVerts(v) = length(Struct.Vdat(v).nverts) == 3;
        end
        bulkVerts = find(goodVerts);

        bulkCells = zeros(1,length(Struct.Cdat));
        extCells = zeros(1,length(Struct.Cdat));
        for c = 1:length(Struct.Cdat)
           if (all(ismember(Struct.Cdat(c).nverts,bulkVerts)))
               bulkCells(c) = 1; 
           elseif (any(ismember(Struct.Cdat(c).nverts,bulkVerts)))
               extCells(c) = 1;
           end
        end

        bulkCells = find(bulkCells);
        extCells = find(extCells);
        involvedCells = [bulkCells,extCells];

        Tri = zeros(length(bulkVerts),3);
        for v = 1:length(bulkVerts)
           nCells = Struct.Vdat(bulkVerts(v)).ncells;
           for c = 1:length(nCells)
               Tri(v,c) = find(involvedCells==nCells(c));
           end
        end

        % Compute vertex/cell positions.
        rV = zeros(length(bulkVerts),2);
        for ii = 1:length(bulkVerts)
            rV(ii,1) = double(Struct.Vdat(bulkVerts(ii)).vertxcoord);
            rV(ii,2) = double(Struct.Vdat(bulkVerts(ii)).vertycoord);
        end

        rC = zeros(length(Struct.Cdat),3);
        for ii = 1:length(Struct.Cdat)
            rC(ii,1:2) = Struct.Cdat(ii).centroid.coord;
        end
        rC = rC(involvedCells,:);

        if (size(Tri,2) == 3)
            [ Tri ] = fitDual.orderTri( rC(:,1:2), Tri );
        end
        
    else
        C = length(Struct.Cdat);
        bulkCells = 1:C;
        
        three_fold = [Struct.Cdat.all_threefold];
        
        if (isfield(Struct.Cdat,'hole'))
            holes = [Struct.Cdat.hole];
        else
            holes = zeros(size(three_fold));
            holes(1) = 1;
        end
        involvedCells = (three_fold==1) & (holes==0);
    
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

        % Bulk verts are verts connected to bulk cells
        bulkVerts = [];
        for bc = bulkCells
            bulkVerts = horzcat(bulkVerts,Struct.Cdat(bc).nverts);
        end
        bulkVerts = unique(bulkVerts);

        % Find extCells
        involvedCells = [];
        for bv = bulkVerts
            involvedCells = horzcat(involvedCells,Struct.Vdat(bv).ncells);
        end

        involvedCells = unique(involvedCells);
        extCells = involvedCells(~ismember(involvedCells,bulkCells));

        involvedCells = [bulkCells,extCells];

        Tri = zeros(length(bulkVerts),3);
        for v = 1:length(bulkVerts)
           nCells = Struct.Vdat(bulkVerts(v)).ncells;
           for c = 1:length(nCells)
               Tri(v,c) = find(involvedCells==nCells(c));
           end
        end
        
        % Compute vertex/cell positions.
        rV = zeros(length(bulkVerts),2);
        for ii = 1:length(bulkVerts)
            rV(ii,1) = double(Struct.Vdat(bulkVerts(ii)).vertxcoord);
            rV(ii,2) = double(Struct.Vdat(bulkVerts(ii)).vertycoord);
        end
        
        rC = zeros(length(Struct.Cdat),3);
        for ii = 1:length(Struct.Cdat)
            rC(ii,1:2) = Struct.Cdat(ii).centroid.coord;
        end
        rC = rC(involvedCells,:);

        [ Tri ] = fitDual.orderTri( rC(:,1:2), Tri );

    end
    
end

