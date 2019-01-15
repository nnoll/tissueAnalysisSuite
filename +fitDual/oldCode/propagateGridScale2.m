function [ scale, goodTiles ] = propagateGridScale2( PN, ERes, L, delta )
    % RESCALE GRID DUAL 
    
    Ec = 10;
%     centerInd = round(size(PN,1)/2) + size(PN,1)*(round(size(PN,2)/2)-1);
%     centerInd = 10 + size(PN,1)*(13-1);
    
    mask = ERes < Ec;
    mask(1:(1+delta),:) = 0;
    mask(:,1:(1+delta)) = 0;
    mask((size(mask,1)-delta):size(mask,1),:) = 0;
    mask(:,(size(mask,2)-delta):size(mask,2)) = 0;

    CC = bwconncomp(mask,4);
    mask = zeros(size(mask));
    mask(CC.PixelIdxList{1}) = 1;
    goodTiles = find(mask);

    nNeighs = mask(goodTiles+1) + mask(goodTiles-1) + mask(goodTiles+size(mask,1)) + mask(goodTiles-size(mask,1));
%     fixedTiles = find(nNeighs == 1);
    
    nTiles = size(PN,1)*size(PN,2);
    delRows = 1:nTiles;
    delRows(goodTiles) = [];

    nEqns = 0;
    qShared = cell(nTiles,nTiles,3);
    
    R0 = regionprops(L,'Area');
    R0 = sqrt([R0.Area]/pi);
    
    nTs = zeros(nTiles,nTiles);
    Zc = 20;
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
                    
                    if (size(q1,1) > Zc)
                        randInd = randsample(1:size(q1,1),Zc);
                        q1 = q1(randInd,:);
                        q2 = q2(randInd,:);
                        ind1 = ind1(randInd);
                    end
                    
                    qShared{indT1,indT2,1} = q1;
                    qShared{indT1,indT2,2} = q2;
                    qShared{indT1,indT2,3} = R0(PN{ii,jj}.cellLabels(ind1))';
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
                    
                    if (size(q1,1) > Zc)
                        randInd = randsample(1:size(q1,1),Zc);
                        q1 = q1(randInd,:);
                        q2 = q2(randInd,:);
                        ind1 = ind1(randInd);
                    end
                    
                    qShared{indT1,indT3,1} = q1;
                    qShared{indT1,indT3,2} = q2;
                    qShared{indT1,indT3,3} = R0(PN{ii,jj}.cellLabels(ind1))';

                    nEqns = nEqns + size(q1,1);
                    nTs(indT1,indT2) = size(q1,1);

                end
            end
        end
    end
    
%     cdfplot(nTs(nTs>0))
%     pause
%     nnTs = sum(nTs,2);
%     nnTs = reshape(nnTs,size(PN,1),size(PN,2));
%     imagesc(nnTs)
%     pause
    
    L = zeros(2*nEqns,3*nTiles);
    b = zeros(2*nEqns,1);
        
    % Fill in matrix elements
    row = 1;
    for n = 1:nTiles
        for nn = 1:nTiles
            if (~isempty(qShared{n,nn}))
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
    
    delInd = delRows;
    delInd = [delInd,delInd+nTiles,delInd+2*nTiles];

    L(:,delInd) = [];
    nT = size(L,2)/3;
    
    cons = zeros(1,size(L,2));
    cons(1,1:nT) = 1/nT;
    cons(2,(nT+1:2*nT)) = 1/nT;
    cons(3,(2*nT+1:3*nT)) = 1/nT;
    
    L = [L;100*cons];
    b = [b;100;0;0];
    
    L = sparse(L);
    Lambda = L \ b;
%     plot(L*Lambda)
%     pause
    
    Lambda = reshape(Lambda,length(Lambda)/3,3);
    
    changedInd = 1:nTiles;
    changedInd(delRows) = [];

    scale = zeros(nTiles,1);
    for n = 1:length(Lambda)
        scale(changedInd(n)) = Lambda(n);
    end
    scale = reshape(scale,size(PN,1),size(PN,2));    
    
end

