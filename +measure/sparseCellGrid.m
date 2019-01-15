function [ xG, yG ] = sparseCellGrid( Struct, nX, nY, d )
    % CELL GRID 

    if (isfield(Struct,'labelMat'))
        L = Struct.labelMat;
        S = regionprops(L ~= 1,'BoundingBox','Area');
        [~,ind] = max([S.Area]);
        XStart = S(ind).BoundingBox(1);
        YStart = S(ind).BoundingBox(2);

        XFinish = S(ind).BoundingBox(1) + S(ind).BoundingBox(3);
        YFinish = S(ind).BoundingBox(2) + S(ind).BoundingBox(4);
    else
        XStart = (min([Struct.Vdat.vertxcoord]));
        XFinish = (max([Struct.Vdat.vertxcoord]));
        YStart = (min([Struct.Vdat.vertycoord]));
        YFinish = (max([Struct.Vdat.vertycoord]));        
    end
    
    XG = linspace(XStart,XFinish,nX+1);
    dX = d*(XG(2)-XG(1));
    xG(1,:) = XG(1:end-1);
    xG(2,:) = XG(2:end);
    
    YG = linspace(YStart,YFinish,nY+1);
    dY = d*(YG(2)-YG(1));
    yG(1,:) = YG(1:end-1);
    yG(2,:) = YG(2:end);
    
    xG(1,2:end) = xG(1,2:end) - dX;
    xG(2,1:(end-1)) = xG(2,1:(end-1)) + dX;

    yG(1,2:end) = yG(1,2:end) - dY;
    yG(2,1:(end-1)) = yG(2,1:(end-1)) + dY;

    if (isfield(Struct,'labelMat'))
       xG = round(xG);
       yG = round(yG);
    end
    
end

