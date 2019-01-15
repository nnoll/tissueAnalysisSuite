function [ xG, yG, cellLabels, plotCell ] = cellGrid( L, nX, nY )
    % CELL GRID 

    S = regionprops(L ~= 1,'BoundingBox','Area');
    [~,ind] = max([S.Area]);
    XStart = S(ind).BoundingBox(1);
    YStart = S(ind).BoundingBox(2);

    XFinish = S(ind).BoundingBox(1) + S(ind).BoundingBox(3);
    YFinish = S(ind).BoundingBox(2) + S(ind).BoundingBox(4);

    XG = linspace(XStart,XFinish,nX+1);
    xCom = .5*(XG(2:end) + XG(1:end-1));
    xG(1,:) = [XG(1:end-1),xCom(1:end-1)];
    xG(2,:) = [XG(2:end),xCom(2:end)];
    [~,ind] = sort(xG(1,:));
    xG(1,:) = xG(1,ind);
    xG(2,:) = xG(2,ind);
    
    YG = linspace(YStart,YFinish,nY+1);
    yCom = .5*(YG(2:end) + YG(1:end-1));
    yG(1,:) = [YG(1:end-1),yCom(1:end-1)];
    yG(2,:) = [YG(2:end),yCom(2:end)];
    [~,ind] = sort(yG(1,:));
    yG(1,:) = yG(1,ind);
    yG(2,:) = yG(2,ind);
    
    if (nargout > 2)
        S = regionprops(L,'Centroid','PixelIdxList');
        rC = vertcat(S.Centroid);
        rC(1,:) = inf;

        cellLabels = cell(nY,nX);
        plotCellR = zeros(size(L));
        plotCellG = zeros(size(L));
        plotCellB = zeros(size(L));

        color = jet(size(xG,2)*size(yG,2));
        color = color(randperm(size(color,1)),:);
        color = reshape(color,size(yG,2),size(xG,2),3);

        for ii = 1:size(yG,2)
            for jj = 1:size(xG,2)
                cellLabels{ii,jj} = find( (rC(:,1) >= xG(1,jj)) & (rC(:,1) <= xG(2,jj)) & ...
                                          (rC(:,2) >= yG(1,ii)) & (rC(:,2) <= yG(2,ii)) );
                plotCellR( vertcat(S(cellLabels{ii,jj}).PixelIdxList) ) = plotCellR( vertcat(S(cellLabels{ii,jj}).PixelIdxList) ) + .25*color(ii,jj,1);
                plotCellG( vertcat(S(cellLabels{ii,jj}).PixelIdxList) ) = plotCellG( vertcat(S(cellLabels{ii,jj}).PixelIdxList) ) + .25*color(ii,jj,2);
                plotCellB( vertcat(S(cellLabels{ii,jj}).PixelIdxList)) = plotCellB( vertcat(S(cellLabels{ii,jj}).PixelIdxList)) + .25*color(ii,jj,3);
            end
        end

        plotCell = cat(3,plotCellR,plotCellG,plotCellB);
    end
    
end

