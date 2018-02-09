function [ scale, goodTiles ] = propagateGridScale( PN, ERes, delta )
    % RESCALE GRID DUAL 
    
    Ec = 10;
%     centerInd = round(size(PN,1)/2) + size(PN,1)*(round(size(PN,2)/2)-1);
%     centerInd = 10 + size(PN,1)*(13-1);
    if (nargin < 3)
        delta = 0;
    end
    
    mask = ERes < Ec;

    if ( delta > 0 )
        delta = delta - 1;
        mask(1:(1+delta),:) = 0;
        mask(:,1:(1+delta)) = 0;
        mask((size(mask,1)-delta):size(mask,1),:) = 0;
        mask(:,(size(mask,2)-delta):size(mask,2)) = 0;
    end

    CC = bwconncomp(mask,4);
    len = cellfun(@length,CC.PixelIdxList);
    [~,ind] = max(len);
    mask = zeros(size(mask));
    mask(CC.PixelIdxList{ind}) = 1;
    
    goodTiles = find(mask);
%     bndryTiles = find(bwperim(mask));
    bndryTiles = [];
%     centerInd = round(size(PN,1)/2) + size(PN,1)*(round(size(PN,2)/2)-1)
    nTiles = size(PN,1)*size(PN,2);
    delRows = 1:nTiles;
    delRows(goodTiles) = [];

    nEqns = 0;
    TShared = cell(nTiles,nTiles,2);

    nTs = zeros(nTiles,nTiles);
    
    for ii = 1:size(PN,1)
        for jj = 1:size(PN,2)
            
            indT1 = ii + size(PN,1)*(jj-1);
            % Get bounds for neighboring tiles.
            yM = min(size(PN,1),ii+1);
            xM = min(size(PN,2),jj+1);
            
            if (ismember(indT1,goodTiles))
                
                indT2 = yM + size(PN,1)*(jj-1);
                indT3 = ii + size(PN,1)*(xM-1);
                
                if (ismember(indT2,goodTiles))
                    Tp1 = sqrt( sum((PN{ii,jj}.d0*PN{ii,jj}.q).^2,2) );
                    Tp2 = sqrt( sum((PN{yM,jj}.d0*PN{yM,jj}.q).^2,2) );

                    ind1 = ismember(PN{ii,jj}.cellLabels,PN{yM,jj}.cellLabels); % All the cells in patch ii,jj that are also in its neighbor.
                    bonds = sum(abs(PN{ii,jj}.d0(:,ind1)),2) == 2; % Record bonds where both cells are still in overlapping geometry.
                    
                    T1 = Tp1(bonds);
                    T2 = zeros(size(T1));
                    
                    n = 1;
                    badBonds = [];
                    for b = find(bonds)'
                       cellInd = (PN{ii,jj}.d0(b,:)~=0);
                       cellLabel = PN{ii,jj}.cellLabels(cellInd);
                       ind2 = ismember(PN{yM,jj}.cellLabels,cellLabel);
                       bond2 = sum(abs(PN{yM,jj}.d0(:,ind2)),2) == 2;
                       if (sum(bond2) == 1)
                           T2(n) = Tp2(bond2);
                       else
                           badBonds = [badBonds,n];
                       end
                       n = n + 1;
                    end
                    T1(badBonds) = [];
                    T2(badBonds) = [];

                    TShared{indT1,indT2,1} = T1;
                    TShared{indT1,indT2,2} = T2;
                    nTs(indT1,indT2) = length(T1);
                    
                    nEqns = nEqns + length(T1);
                end

                if (ismember(indT3,goodTiles))

                    Tp1 = sqrt( sum((PN{ii,jj}.d0*PN{ii,jj}.q).^2,2) );
                    Tp2 = sqrt( sum((PN{ii,xM}.d0*PN{ii,xM}.q).^2,2) );

                    ind1 = ismember(PN{ii,jj}.cellLabels,PN{ii,xM}.cellLabels); % All the cells in patch ii,jj that are also in its neighbor.
                    bonds = sum(abs(PN{ii,jj}.d0(:,ind1)),2) == 2; % Record bonds where both cells are still in overlapping geometry.
                    
                    T1 = Tp1(bonds);
                    T2 = zeros(size(T1));
                    
                    n = 1;
                    badBonds = [];
                    for b = find(bonds)'
                       cellInd = (PN{ii,jj}.d0(b,:)~=0);
                       cellLabel = PN{ii,jj}.cellLabels(cellInd);
                       ind2 = ismember(PN{ii,xM}.cellLabels,cellLabel);
                       bond2 = sum(abs(PN{ii,xM}.d0(:,ind2)),2) == 2;
                       if (sum(bond2) == 1)
                           T2(n) = Tp2(bond2);
                       else
                           badBonds = [badBonds,n]; 
                       end
                       n = n + 1;
                    end
                    T1(badBonds) = [];
                    T2(badBonds) = [];

                    TShared{indT1,indT3,1} = T1;
                    TShared{indT1,indT3,2} = T2;
                    nTs(indT1,indT3) = length(T1);

                    nEqns = nEqns + length(T1);

                end
            end
        end
    end

    L = zeros(nEqns,nTiles);
    b = zeros(nEqns,1);

    % Fill in matrix elements
    row = 1;
    for n = 1:nTiles
        for nn = 1:nTiles
            if (~isempty(TShared{n,nn}))
                if (ismember(n,bndryTiles) && ~ismember(nn,bndryTiles))
                    rowEnd = row + size(TShared{n,nn,1},1) - 1;
                    L(row:rowEnd,nn) = TShared{n,nn,2};
                    b(row:rowEnd) = TShared{n,nn,1};
                    
                    row = rowEnd + 1;
                elseif (ismember(nn,bndryTiles) && ~ismember(n,bndryTiles))
                    rowEnd = row + size(TShared{n,nn,1},1) - 1;
                    L(row:rowEnd,n) = TShared{n,nn,1};
                    b(row:rowEnd) = TShared{n,nn,2};

                    row = rowEnd + 1;
                else
                    rowEnd = row + size(TShared{n,nn,1},1) - 1;
                    L(row:rowEnd,n) = TShared{n,nn,1};
                    L(row:rowEnd,nn) = -TShared{n,nn,2};

                    row = rowEnd + 1;
                end
            end
        end
    end

    delInd = [delRows,bndryTiles'];

    L(:,delInd) = [];
    L = [L;ones(1,size(L,2))];
    b = [b;size(L,2)];
    
    L = sparse(L);
    Lambda = L \ b;     

    changedInd = 1:nTiles;
    changedInd(delInd) = [];

    scale = zeros(nTiles,1);
    for n = 1:length(Lambda)
        scale(changedInd(n)) = Lambda(n);
    end
    scale(bndryTiles') = 1;
    scale = reshape(scale,size(PN,1),size(PN,2));
    
end

