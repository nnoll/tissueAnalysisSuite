function [ dC, dV, cAdj, involvedVerts, involvedCells ] = computeSubDiffOperators( Struct, XRange, YRange )
    %COMPUTE DIFF OPERATORS 

    [ ~, bulkCells, extCells, bulkVerts ] = fitDual.returnSubGraph( Struct, XRange, YRange );
    [ extVerts ] = fitDual.returnExtVerts( Struct, bulkVerts );

    involvedVerts = [bulkVerts,extVerts];
    involvedCells = [bulkCells,extCells];
    
%     plot.skel(Struct,'r',1)
%     hold on
%     scatter([Struct.Vdat(bulkVerts).vertxcoord],[Struct.Vdat(bulkVerts).vertycoord],'b','filled');
%     pause

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
%     plot.skel(Struct,'b',0);
%     hold all
%     scatter([Struct.Vdat(bulkVerts).vertxcoord],[Struct.Vdat(bulkVerts).vertycoord],'g','filled')
%     scatter([Struct.Vdat(extVerts).vertxcoord],[Struct.Vdat(extVerts).vertycoord],'r','filled')
    
    for ii = 1:nCells
        neighCells = find(cAdj(ii,:));
        neighCells = neighCells(neighCells>ii);
        for nc = neighCells
           dC(nB,ii) = 1;
           dC(nB,nc) = -1;

           c1 = involvedCells(ii);
           c2 = involvedCells(nc);
           
           bVerts = Struct.Cdat(c1).nverts(ismember(Struct.Cdat(c1).nverts,Struct.Cdat(c2).nverts));

           if (length(bVerts) == 2 && sum(ismember(bVerts,involvedVerts)) == 2)
               dV(nB,involvedVerts==bVerts(1)) = 1;
               dV(nB,involvedVerts==bVerts(2)) = -1;
%                scatter([Struct.Vdat(bVerts).vertxcoord],[Struct.Vdat(bVerts).vertycoord],'ko','LineWidth',2)

%                if (sum(ismember(bVerts,extVerts)) >= 1)
%                    clf
%                    plot.skel(Struct,'b',0);
%                    rc = [Struct.Cdat.centroid];
%                    rc = vertcat(rc.coord);
%                    hold all
%                    nB
%                    scatter([Struct.Vdat(bulkVerts).vertxcoord],[Struct.Vdat(bulkVerts).vertycoord],'go','LineWidth',2)
%                    scatter([Struct.Vdat(extVerts).vertxcoord],[Struct.Vdat(extVerts).vertycoord],'ko','LineWidth',2)
%                    scatter([Struct.Vdat(bVerts).vertxcoord],[Struct.Vdat(bVerts).vertycoord],'r','filled')
%                    scatter(rc([c1,c2],1),rc([c1,c2],2),'c','filled')
%                    set(gca,'XLim',[350,650],'YLim',[250,550])
%                    pause(.1)
%                end
           end
           
           nB = nB + 1;
           
        end
    end
%         pause
    badVerts = sum(abs(dV),1) == 0;
    dV(:,badVerts) = [];
    involvedVerts(badVerts) = [];

%     figure(2)
    badBonds = sum(abs(dV),2) == 0;
    dC(badBonds,:) = [];
    dV(badBonds,:) = [];
%     plot(sum(abs(dV),2))
%     hold all
%     plot(sum(abs(dC),2))
%     find(sum(abs(dV(:,end-length(extVerts):end)),2)==2)
%     pause
end

