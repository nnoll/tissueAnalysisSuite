function [ PN, goodTiles ] = rescaleGridDual( PN, ERes, alpha )
    % RESCALE GRID DUAL 
    
    if (nargin == 2)
        alpha = 10;
    end
    
    Ec = 10;
    mask = ERes < Ec;
    mask(1,:) = 0;
    mask(:,1) = 0;
    mask(size(mask,1),:) = 0;
    mask(:,size(mask,2)) = 0;
   
%     mask = imerode(mask,strel('disk',1));

    CC = bwconncomp(mask,4);
    mask = zeros(size(mask));
    mask(CC.PixelIdxList{1}) = 1;

    goodTiles = find(mask);
    
    nTiles = size(PN,1)*size(PN,2);
    delRows = 1:nTiles;
    delRows(goodTiles) = [];
    
    qShared = cell(nTiles,nTiles,2);
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

                    ind1 = ismember(PN{ii,jj}.cellLabels,PN{yM,jj}.cellLabels);
                    q1 = PN{ii,jj}.q(ind1,:);
                    q2 = zeros(size(q1));
                    
                    n = 1;
                    for ind2 = find(ind1)
                       q2(n,:) = PN{yM,jj}.q(ismember(PN{yM,jj}.cellLabels,PN{ii,jj}.cellLabels(ind2)),:);
                       n = n + 1;
                    end
                    
%                     ii
%                     jj
%                     yM
%                     figure(1)
%                     clf
%                     PN{ii,jj}.plotPrimal('b')
%                     hold all
%                     PN{yM,jj}.plotPrimal('r')
%                     scatter(PN{ii,jj}.q(:,1),PN{ii,jj}.q(:,2),'bo')
%                     scatter(q1(:,1),q1(:,2),'c','filled')
%                     scatter(PN{yM,jj}.q(:,1),PN{yM,jj}.q(:,2),'ro')
%                     scatter(q2(:,1),q2(:,2),'y','filled')
%                     quiver(q1(:,1),q1(:,2),q2(:,1)-q1(:,1),q2(:,2)-q1(:,2),0,'k')
%                     
%                     figure(2)
%                     scatter(q1(:,1),q2(:,1))
%                     
%                     figure(3)
%                     scatter(q1(:,2),q2(:,2))
%                     pause
                    
                    qShared{indT1,indT2,1} = q1;
                    qShared{indT1,indT2,2} = q2;
                    nEqns = nEqns + size(q1,1);
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

                end
            end
        end
    end
    
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
    
    delInd = [delRows,delRows+nTiles,delRows+2*nTiles];
    L(:,delInd) = [];
    
    nTs = size(L,2)/3;
%     plot( L* [ones(nTs,1);zeros(nTs,1);zeros(nTs,1)] )
%     pause
    
    L = [L; [alpha*median(mean(abs(L),2))*ones(1,nTs),zeros(1,nTs),zeros(1,nTs)]; ...
        [zeros(1,nTs),ones(1,nTs),zeros(1,nTs)]; ...
        [zeros(1,nTs),zeros(1,nTs),ones(1,nTs)]];
    b = [b;alpha*median(mean(abs(L),2))*nTs;0;0];

    L = sparse(L);
%     Lambda = L \ b;
%     A = -[eye(nTs),zeros(nTs),zeros(nTs)];
%     d = zeros(size(A,1),1);
    lb = [-10*ones(nTs,1);-inf*ones(nTs,1);-inf*ones(nTs,1)];
    ub = [10*ones(nTs,1);inf*ones(nTs,1);inf*ones(nTs,1)];
    Lambda = lsqlin(L,b,[],[],[],[],lb,ub);
    Lambda = reshape(Lambda,length(Lambda)/3,3);
    
    plot(Lambda(:,1))
    changedInd = 1:nTiles;
    changedInd(delRows) = [];

    n = 1;
    for ii = changedInd
%         if (n == 329)
%             ii
%             [x,y] = ind2sub([size(PN,1),size(PN,2)],ii)
%             Lambda(n,1)
%         end
        if (~isempty(PN{ii}) && Lambda(n,1) > 0)
            PN{ii} = PN{ii}.rescaleDual(Lambda(n,1),[0,0]); 
        end
        n = n + 1;
    end
    
end

