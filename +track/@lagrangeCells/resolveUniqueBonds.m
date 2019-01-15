function [ bonds ] = resolveUniqueBonds( this, Struct )
    % RETURN UNIQUE BONDS 

    bonds(1) = struct('label',[],'cells',[],'t0',[]);
    
    % Initialize bond data structure at t=1
    nb = 1;
    [ IDs, labels ] = this.returnIDMap( 1 );
    for b = 1:length(Struct(1).Bdat)
        if (all(ismember(Struct(1).Bdat(b).cells,labels)))
            bonds(nb).label(1) = b;
            bonds(nb).cells = [IDs(labels==Struct(1).Bdat(b).cells(1)),IDs(labels==Struct(1).Bdat(b).cells(2))];
            bonds(nb).t0 = 1;
            nb = nb + 1;
        end
    end
    
    for t = 2:length(Struct)
        [ IDs, labels ] = this.returnIDMap( t );
        ind = zeros(length(Struct(t).Bdat),1);
        ind(labels) = IDs;
        bCells = vertcat(Struct(t).Bdat.cells); 
        bCells = ind(bCells);
        
        uBCells = vertcat(bonds.cells);
        [~,ind,indU] = intersect(bCells,uBCells,'rows');
        for ii = 1:length(indU)
           bonds(indU(ii)).label(t) = ind(ii);
        end
        
        % Add any bonds that were not matched and have both IDs
        % represented.
        newBonds = 1:size(bCells,1);
        newBonds = newBonds(~ismember(newBonds,ind));
        newBonds = newBonds(sum(bCells(newBonds,:)~=0,2) == 2);
        n = length(bonds) + 1;
        for ii = 1:length(newBonds)
            bonds(n).label(t) = newBonds(ii);
            bonds(n).cells = bCells(newBonds(ii),:);
            bonds(n).t0 = t;
            n = n + 1;
        end
    end
    
end

