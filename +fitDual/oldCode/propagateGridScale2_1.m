function [ scale, goodTiles ] = propagateGridScale2( PN, ERes, delta )
    % RESCALE GRID DUAL 
    
    Ec = 10;
    if (nargin < 3)
        delta = 0;
    end
    
    mask = ERes < Ec;
    if ( delta > 0 )
        mask(1:delta,:) = 0;
        mask(:,1:delta) = 0;
        mask((size(mask,1)-delta+1):size(mask,1),:) = 0;
        mask(:,(size(mask,2)-delta+1):size(mask,2)) = 0;
    end
    
    CC = bwconncomp(mask,4);
    len = cellfun(@length,CC.PixelIdxList);
    [~,ind] = max(len);

    mask = zeros(size(mask));
    mask(CC.PixelIdxList{ind}) = 1;
    goodTiles = find(mask);

    nTiles = size(PN,1)*size(PN,2);
    delRows = 1:nTiles;
    delRows(goodTiles) = [];

    nEqns = 0;
    qShared = cell(nTiles,nTiles,2);
 
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
                rowEnd = row + size(qShared{n,nn,1},1) - 1;
                Qscale = sqrt( sqrt( sum(qShared{n,nn,1}.^2,2) .* sum(qShared{n,nn,2}.^2,2) ) );
%                 Qscale = ones(size(Qscale));

                % x component
                L(row:rowEnd,n) = qShared{n,nn,1}(:,1)./Qscale;
                L(row:rowEnd,n+nTiles) = 1./Qscale;
                L(row:rowEnd,nn) = -qShared{n,nn,2}(:,1)./Qscale;
                L(row:rowEnd,nn+nTiles) = -1./Qscale;
                b(row:rowEnd) = (-qShared{n,nn,1}(:,1) + qShared{n,nn,2}(:,1))./Qscale;

                % y component
                L(nEqns+(row:rowEnd),n) = qShared{n,nn,1}(:,2)./Qscale;
                L(nEqns+(row:rowEnd),n+2*nTiles) = 1./Qscale;
                L(nEqns+(row:rowEnd),nn) = -qShared{n,nn,2}(:,2)./Qscale;
                L(nEqns+(row:rowEnd),nn+2*nTiles) = -1./Qscale;
                b(nEqns+(row:rowEnd)) = (-qShared{n,nn,1}(:,2) + qShared{n,nn,2}(:,2))./Qscale;

                row = rowEnd + 1;
            end
        end
    end
    
    delInd = delRows;
    delInd = [delInd,delInd+nTiles,delInd+2*nTiles];

    L(:,delInd) = [];
    nT = size(L,2)/3;
    cons = zeros(3,size(L,2));

    cons(1,(1):1*nT) = 1/nT;
    cons(2,(nT+1):2*nT) = 1/nT;
    cons(3,(2*nT+1):3*nT) = 1/nT;

    L = [L;10*cons];
    b = [b;0;0;0];

    L = sparse(L);
    Lambda = L \ b;
    
    Lambda = reshape(Lambda,length(Lambda)/3,3);
    Lambda = Lambda(:,1);

    changedInd = 1:nTiles;
    changedInd(delRows) = [];

    scale = zeros(nTiles,1);
    for n = 1:length(Lambda)
        scale(changedInd(n)) = Lambda(n);
    end

    scale = reshape(scale,size(PN,1),size(PN,2));    
    scale = abs(1 + scale);

end

