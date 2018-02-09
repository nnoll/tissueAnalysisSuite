function [ Tri, bulkCells, extCells, bulkVerts, rC, rV ] = returnSubGraph( Struct, XRange, YRange )
    % FIND SUB GRAPH
    
    C = length(Struct.Cdat);
    bulkCells = 1:C;
    three_fold = [Struct.Cdat.all_threefold];
    if (isfield(Struct.Cdat,'hole'))
        holes = [Struct.Cdat.hole];
    else
        holes = zeros(size(bulkCells));
        holes(1) =1;
    end
%     Rc = [Struct.Cdat.centroid];
%     Rc = vertcat(Rc.coord);

    if (isfield(Struct,'labelMat'))
        L = Struct.labelMat;
        containedCells = L(YRange(1):YRange(2),XRange(1):XRange(2));
        containedCells = unique(containedCells(containedCells>1));
    else
        rV = horzcat([Struct.Vdat.vertxcoord]',[Struct.Vdat.vertycoord]');
        containedVerts = ( rV(:,1) > XRange(1) ) & ( rV(:,1) < XRange(2) ) & ...
                         ( rV(:,2) > YRange(1) ) & ( rV(:,2) < YRange(2) );

        containedCells = [];
        for v = find(containedVerts)
            containedCells = [containedCells,Struct.Vdat(v).ncells];
        end
        containedCells = unique(containedCells);
    end

    involvedCells = (three_fold==1) & (holes==0) & ismember(bulkCells,containedCells);

    % Find bndryCells
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
    
%     plot.skel(Struct,'r',0)
%     hold on
%     scatter(rV(:,1),rV(:,2),'b','filled')
%     scatter(rC(:,1),rC(:,2),'g','filled')
%     pause
%     plot.skel(Struct,'r',1)
%     hold on
%     scatter(rV(:,1),rV(:,2),'b','filled');
%     scatter(rC(:,1),rC(:,2),'g','filled');
%     rectangle('Position',[XRange(1),YRange(1),XRange(2)-XRange(1),YRange(2)-YRange(1)],'EdgeColor','k','FaceColor','none','LineWidth',2)
%     pause

    [ Tri ] = fitDual.orderTri( rC(:,1:2), Tri );

end

