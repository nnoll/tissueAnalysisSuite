function [ PN ] = rescaleGridDual( PN, ERes )
    % RESCALE GRID DUAL 
    
    Ec = 10;
    centerInd = round(size(PN,1)/2) + size(PN,1)*(round(size(PN,2)/2)-1);
    nTiles = size(PN,1)*size(PN,2);
    
    qShared = cell(nTiles,nTiles,2);
    delRows = [];
    nEqns = 0;
    for ii = 1:size(PN,1)
        for jj = 1:size(PN,2)
            indT1 = ii + size(PN,1)*(jj-1);
            % Get bounds for neighboring tiles.
            yM = min(size(PN,1),ii+1);
            xM = min(size(PN,2),jj+1);
            
            if (isempty(PN{ii,jj}) || ERes(ii,jj) >= Ec)
                delRows = [delRows,indT1];
            end
            
            if (~isempty(PN{ii,jj}) && ERes(ii,jj) < Ec)

                if (yM > ii && ~isempty(PN{yM,jj}) && ERes(yM,jj) < Ec)

                    ind1 = ismember(PN{ii,jj}.cellLabels,PN{yM,jj}.cellLabels);
                    indT2 = yM + size(PN,1)*(jj-1);

                    q1 = PN{ii,jj}.q(ind1,:);
                    q2 = zeros(size(q1));
                    n = 1;
                    for ind2 = find(ind1)
                       q2(n,:) = PN{yM,jj}.q(ismember(PN{yM,jj}.cellLabels,PN{ii,jj}.cellLabels(ind2)),:);
                       n = n + 1;
                    end
                    
%                     clf
%                     PN{ii,jj}.plotPrimal('c')
%                     hold all
%                     scatter(q1(:,1),q1(:,2),'bo')
%                     PN{yM,jj}.plotPrimal('g')
%                     scatter(q2(:,1),q2(:,2),'r','filled')
%                     pause
                    
                    qShared{indT1,indT2,1} = q1;
                    qShared{indT1,indT2,2} = q2;
                    nEqns = nEqns + size(q1,1);
                end

                if (xM > jj && ~isempty(PN{ii,xM}) && ERes(ii,xM) < Ec)

                    ind1 = ismember(PN{ii,jj}.cellLabels,PN{ii,xM}.cellLabels);
                    indT2 = ii + size(PN,1)*(xM-1);

                    q1 = PN{ii,jj}.q(ind1,:);
                    q2 = zeros(size(q1));
                    n = 1;
                    for ind2 = find(ind1)
                       q2(n,:) = PN{ii,xM}.q(ismember(PN{ii,xM}.cellLabels,PN{ii,jj}.cellLabels(ind2)),:);
                       n = n + 1;
                    end
                    
                    qShared{indT1,indT2,1} = q1;
                    qShared{indT1,indT2,2} = q2;
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
                if (n == centerInd)
                    rowEnd = row + size(qShared{n,nn,1},1) - 1;
                    
                    % x component
                    L(row:rowEnd,nn) = qShared{n,nn,2}(:,1);
                    L(row:rowEnd,nn+nTiles) = 1;
                    b(row:rowEnd) = qShared{n,nn,1}(:,1);
                    
                    % y component
                    L(nEqns+(row:rowEnd),nn) = qShared{n,nn,2}(:,2);
                    L(row:rowEnd,nn+2*nTiles) = 1;
                    b(nEqns+(row:rowEnd)) = qShared{n,nn,1}(:,2);
                    row = rowEnd + 1;
                    
                elseif (nn == centerInd)
                    rowEnd = row + size(qShared{n,nn,1},1) - 1;
                    
                    % x component
                    L(row:rowEnd,n) = qShared{n,nn,1}(:,1);
                    L(row:rowEnd,n+nTiles) = 1;
                    b(row:rowEnd) = qShared{n,nn,2}(:,1);
                    
                    % y component
                    L(nEqns+(row:rowEnd),n) = qShared{n,nn,1}(:,2);
                    L(row:rowEnd,n+2*nTiles) = 1;
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
    delInd = [delRows,centerInd];
    delInd = [delInd,delInd+nTiles,delInd+2*nTiles];
    plot(sum(abs(L),1))
    pause
    L(:,delInd) = [];
    plot(sum(abs(L),1))
    pause
    N = null(L);
    size(N)
    for ii = 1:size(N,2)
        plot(N(:,ii))
        pause
    end
    L = sparse(L);
    Lambda = L \ b;
    Lambda = reshape(Lambda,length(Lambda)/3,3);
    
    changedInd = 1:nTiles;
    changedInd([centerInd,delRows]) = [];

    n = 1;
    for ii = changedInd
        if (~isempty(PN{ii}) && Lambda(n,1) > 0)
            PN{ii} = PN{ii}.rescaleDual(Lambda(n,1),[0,0]); 
        end
        n = n + 1;
    end
    
end

