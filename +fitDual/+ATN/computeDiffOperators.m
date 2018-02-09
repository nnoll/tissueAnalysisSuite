function [ dC, dV, cAdj, involvedVerts, involvedCells ] = computeDiffOperators( Struct, extCell )
    %COMPUTE DIFF OPERATORS 

    if (nargin == 1)
        extCell = 0;
    end
    
    [ ~, bulkCells, extCells, bulkVerts ] = fitDual.returnGraph( Struct, extCell );
    [ extVerts ] = fitDual.returnExtVerts( Struct, bulkVerts );
        
    involvedVerts = [bulkVerts,extVerts];
    involvedCells = [bulkCells,extCells];
    
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
    
    badBonds = sum(abs(dV),2) < 2;
    
    dC(badBonds,:) = [];
    dV(badBonds,:) = [];

end

