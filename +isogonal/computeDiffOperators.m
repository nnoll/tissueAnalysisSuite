function [ dC, dV, cAdj, involvedVerts, involvedCells ] = computeDiffOperators( Struct )
    %COMPUTE DIFF OPERATORS 

    C = length(Struct.Cdat);
    
    % Remove external cell and its neighbors
    bulkCells = 1:C;
    bulkCells(bulkCells == 1) = [];
    bulkCells(ismember(bulkCells,Struct.Cdat(1).ncells)) = [];
    bad_cells = [];
    for ii=1:length(bulkCells)
        if (sum(ismember(Struct.Cdat(bulkCells(ii)).ncells,bulkCells))==0)
            bad_cells = horzcat(bad_cells,ii);
        end
    end
    bulkCells(bad_cells) = [];
    
    %Find involved vertices 
    bulkVerts = [];
    for bc = bulkCells
        bulkVerts = horzcat(bulkVerts,Struct.Cdat(bc).nverts);
    end
    
    [bulkVerts,~] = isogonal.count_unique(bulkVerts);
    
    involvedCells = [];
    involvedVerts = [];
    for bv = bulkVerts
        involvedCells = horzcat(involvedCells,Struct.Vdat(bv).ncells);
        involvedVerts = horzcat(involvedVerts,Struct.Vdat(bv).nverts);
    end
    involvedCells = fitDual.count_unique(involvedCells);
    involvedVerts = fitDual.count_unique(involvedVerts);
    
    extCells = involvedCells(~ismember(involvedCells,bulkCells))';
    extVerts = involvedVerts(~ismember(involvedVerts,bulkVerts))';
    
    involvedCells = [bulkCells,extCells];
    involvedVerts = [bulkVerts',extVerts];
    
    nVerts = length(involvedVerts);
    nCells = length(involvedCells);

    %Build cell adjacency matrix.
    cAdj = zeros(nCells);
    nBonds = 0;
    for ii = 1:nCells
        c = involvedCells(ii);
        neighCells = Struct.Cdat(c).ncells;
        for nc = neighCells
            index = find(involvedCells == nc,1);
            if ( ~isempty(index) && cAdj(ii,index) == 0 )
                cAdj(ii,index) = 1;
                cAdj(index,ii) = 1;
                nBonds = nBonds + 1;
            end
        end
    end
    
    % Build difference operators.
    dC = zeros(nBonds,nCells);
    dV = zeros(nBonds,nVerts);
    
    nB = 1;
    for ii = 1:nCells
        neighCells = find(cAdj(ii,:));
        neighCells = neighCells(neighCells>ii);
        for nc = neighCells
           dC(nB,ii) = 1;
           dC(nB,nc) = -1;
           
           c1 = involvedCells(ii);
           c2 = involvedCells(nc);
           
           bVerts = Struct.Cdat(c1).nverts(ismember(Struct.Cdat(c1).nverts,Struct.Cdat(c2).nverts));
           if (length(bVerts) == 2)
               dV(nB,involvedVerts==bVerts(1)) = 1;
               dV(nB,involvedVerts==bVerts(2)) = -1;
           end
           
           nB = nB + 1;
           
        end
    end

    badVerts = sum(abs(dV),1) == 0;
    dV(:,badVerts) = [];
    involvedVerts(badVerts) = [];
    
    badBonds = sum(abs(dV),2) == 0;
    dC(badBonds,:) = [];
    dV(badBonds,:) = [];
    
    % Remove bad cells
    badCells = sum(abs(dC),1) <= 2;
    badBonds = sum(abs(dC(:,badCells)),2) >= 1;
    
    dC(:,badCells) = [];
    dC(badBonds,:) = [];
    dV(badBonds,:) = [];
    
    involvedCells(badCells) = [];
    cAdj(badCells,:) = [];
    cAdj(:,badCells) = [];
    
    badVerts = sum(abs(dV),1) == 0;
    dV(:,badVerts) = [];
    involvedVerts(badVerts) = [];
    
end

