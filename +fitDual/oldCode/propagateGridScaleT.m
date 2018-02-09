function [ scale, goodTiles ] = propagateGridScale( PN, ERes, delta )
    % RESCALE GRID DUAL 
    
    if (nargin == 2)
        delta = 0;
    end
    
    Ec = 10;
    mask = ERes < Ec;
    mask(1:(1+delta),:) = 0;
    mask(:,1:(1+delta)) = 0;
    mask((size(mask,1)-delta):size(mask,1),:) = 0;
    mask(:,(size(mask,2)-delta):size(mask,2)) = 0;

%     mask = imerode(mask,strel('disk',1));

    CC = bwconncomp(mask,4);
    mask = zeros(size(mask));
    mask(CC.PixelIdxList{1}) = 1;
%     mask(:,1:8) = 0;
    
    goodTiles = find(mask);
    
    % Find any tiles that are only connected to one other tile. Fix their
    % scale to be the same as the one its connected to.
    nNeighs = mask(goodTiles+1) + mask(goodTiles-1) + mask(goodTiles+size(mask,1)) + mask(goodTiles-size(mask,1));
%     rgb(:,:,1) = mask;
%     tmp = zeros(size(mask));
%     tmp(goodTiles(nNeighs==1)) = 1;
%     rgb(:,:,2) = tmp;
%     rgb(:,:,3) = 0;
   
    fixedTiles = find(nNeighs == 1);
%     for ii = 1:fixedTiles
%         if (mask(goodTiles(fixedTiles(f))+1) == 1)
%             pairedTiles(ii) = find(goodTiles == goodTiles(fixedTiles(f))+1);
%         elseif (mask(goodTiles(fixedTiles(f))-1) == 1)
%             pairedTiles(ii) = find(goodTiles == goodTiles(fixedTiles(f))-1);
%         elseif (mask(goodTiles(fixedTiles(f))+size(mask,1)) == 1)
%             pairedTiles(ii) = find(goodTiles == goodTiles(fixedTiles(f))+size(mask,1));
%         else
%             pairedTiles(ii) = find(goodTiles == goodTiles(fixedTiles(f))-size(mask,1));
%         end
%     end
    
    nTiles = size(PN,1)*size(PN,2);
    delRows = 1:nTiles;
    delRows(goodTiles) = [];

    TShared = cell(nTiles,nTiles,2);
    nEqns = 0;
    
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

%                     badBonds = ( (T1 < 1) | (T2 < 1) );
%                     T1(badBonds) = [];
%                     T2(badBonds) = [];

                    TShared{indT1,indT2,1} = T1;
                    TShared{indT1,indT2,2} = T2;
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
%                     badBonds = ( (T1 < 1) | (T2 < 1) );
%                     T1(badBonds) = [];
%                     T2(badBonds) = [];

                    TShared{indT1,indT3,1} = T1;
                    TShared{indT1,indT3,2} = T2;
                    
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
                rowEnd = row + size(TShared{n,nn,1},1) - 1;

                L(row:rowEnd,n) = TShared{n,nn,1};
                L(row:rowEnd,nn) = -TShared{n,nn,2};

                row = rowEnd + 1;
            end
        end
    end

    L(:,delRows) = [];
    nTs = size(L,2);
      
%     ASCALE = median(sum(abs(L),2));
    ASCALE = .001;
    L = [L; ASCALE*ones(1,nTs)/nTs];
    b = [b;ASCALE*1];
    
    constraint = zeros(length(fixedTiles),size(L,2));
    for ii = 1:length(fixedTiles)
       constraint(ii,fixedTiles(ii)) = 1; 
    end
    
    L = [L;ASCALE*constraint];
    b = [b;ASCALE*ones(size(constraint,1),1)];
    
    L = sparse(L);
    Lambda = L \ b;
    
    changedInd = 1:nTiles;
    changedInd(delRows) = [];

    scale = zeros(nTiles,1);
    for n = 1:length(Lambda)
        scale(changedInd(n)) = Lambda(n);
    end
    scale = reshape(scale,size(PN,1),size(PN,2));
    
end

