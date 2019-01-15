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
    bndryTiles = find(bwperim(mask));

    nTiles = size(PN,1)*size(PN,2);
    delRows = 1:nTiles;
    delRows(goodTiles) = [];

    nEqns = 0;
    qShared = cell(nTiles,nTiles,3);

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

                    ind1 = ismember(PN{ii,jj}.cellLabels,PN{yM,jj}.cellLabels);
                    q1 = PN{ii,jj}.q(ind1,:);
                    q2 = zeros(size(q1));
                    n = 1;
                    for ind2 = find(ind1)
                       q2(n,:) = PN{yM,jj}.q(ismember(PN{yM,jj}.cellLabels,PN{ii,jj}.cellLabels(ind2)),:);
                       n = n + 1;
                    end
                    
                    qShared{indT1,indT2,1} = q1;
                    qShared{indT1,indT2,2} = q2;
                    
                    nEqns = nEqns + size(q1,1);
                    nTs(indT1,indT2) = size(q1,1);
                end

                if (ismember(indT3,goodTiles))

                    ind1 = ismember(PN{ii,jj}.cellLabels,PN{ii,xM}.cellLabels);

                    q1 = PN{ii,jj}.q(ind1,:);
                    q2 = zeros(size(q1));
                    n = 1;
                    for ind2 = find(ind1)
                       q2(n,:) = PN{ii,xM}.q(ismember(PN{ii,xM}.cellLabels,PN{ii,jj}.cellLabels(ind2)),:);
                       n = n + 1;
                    end
                    
                    qShared{indT1,indT3,1} = q1;
                    qShared{indT1,indT3,2} = q2;

                    nEqns = nEqns + size(q1,1);
                    nTs(indT1,indT2) = size(q1,1);

                end
            end
        end
    end
    
    L = zeros(2*nEqns,3*nTiles);
    b = zeros(2*nEqns,1);
        
    row = 1;
    for n = 1:nTiles
        for nn = 1:nTiles
            if (~isempty(qShared{n,nn}))
                if (ismember(n,bndryTiles) && ~ismember(nn,bndryTiles))
                    rowEnd = row + size(qShared{n,nn,1},1) - 1;
                    
                    % x component
                    L(row:rowEnd,nn) = qShared{n,nn,2}(:,1);
                    L(row:rowEnd,nn+nTiles) = 1;
                    b(row:rowEnd) = qShared{n,nn,1}(:,1);
                    
                    % y component
                    L(nEqns+(row:rowEnd),nn) = qShared{n,nn,2}(:,2);
                    L(nEqns+(row:rowEnd),nn+2*nTiles) = 1;
                    b(nEqns+(row:rowEnd)) = qShared{n,nn,1}(:,2);
                    row = rowEnd + 1;
                    
                elseif (ismember(nn,bndryTiles) && ~ismember(n,bndryTiles))
                    rowEnd = row + size(qShared{n,nn,1},1) - 1;
                    
                    % x component
                    L(row:rowEnd,n) = qShared{n,nn,1}(:,1);
                    L(row:rowEnd,n+nTiles) = 1;
                    b(row:rowEnd) = qShared{n,nn,2}(:,1);
                    
                    % y component
                    L(nEqns+(row:rowEnd),n) = qShared{n,nn,1}(:,2);
                    L(nEqns+(row:rowEnd),n+2*nTiles) = 1;
                    b(nEqns+(row:rowEnd)) = qShared{n,nn,2}(:,2);
                    row = rowEnd + 1;
                else
                    rowEnd = row + size(qShared{n,nn,1},1) - 1;
                    
                    % x component
                    L(row:rowEnd,n) = qShared{n,nn,1}(:,1);
                    L(row:rowEnd,n+nTiles) = 1;
                    L(row:rowEnd,nn) = -qShared{n,nn,2}(:,1);
                    L(row:rowEnd,nn+nTiles) = -1;
                    
                    % y component
                    L(nEqns+(row:rowEnd),n) = qShared{n,nn,1}(:,2);
                    L(nEqns+(row:rowEnd),n+2*nTiles) = 1;
                    L(nEqns+(row:rowEnd),nn) = -qShared{n,nn,2}(:,2);
                    L(nEqns+(row:rowEnd),nn+2*nTiles) = -1;
                    
                    row = rowEnd + 1;
                end
            end
        end
    end
    
    delInd = [delRows,bndryTiles'];
    delInd = [delInd,delInd+nTiles,delInd+2*nTiles];

    L(:,delInd) = [];
    L = sparse(L);
    Lambda = L \ b;
    
    Lambda = reshape(Lambda,length(Lambda)/3,3);
    Lambda = Lambda(:,1);

    changedInd = 1:nTiles;
    changedInd([delRows,bndryTiles']) = [];

    scale = zeros(nTiles,1);
    for n = 1:length(Lambda)
        scale(changedInd(n)) = Lambda(n);
    end

    scale(bndryTiles') = 1;
    scale = abs(reshape(scale,size(PN,1),size(PN,2)));    
    
end

