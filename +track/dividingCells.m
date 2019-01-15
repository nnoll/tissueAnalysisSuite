function [ mitCells ] = dividingCells( L, divCells )
    % DIVIDING CELLS 

    nC = 1;
    for t = 2:size(divCells,4)
        Lt = L(:,:,t-1);
        mitPix = divCells(:,:,3,t) - divCells(:,:,1,t);
        S = regionprops(Lt,'Area','PixelIdxList');
%         rgb(:,:,1) = imdilate(Lt == 0,strel('disk',2));
%         rgb(:,:,2) = mitPix;
%         tmp = zeros(size(rgb,1),size(rgb,2));
        
        implicatedCells = Lt(mitPix==255);
        implicatedCells = implicatedCells(implicatedCells>1);
        uCells = unique(implicatedCells);
        for c = 1:length(uCells)
           nPix = sum(implicatedCells==uCells(c));
           if (nPix >= .75*S(uCells(c)).Area)
               mitCells(nC).label = track.oneCell( L, uCells(c), t-1 );
               mitCells(nC).t0 = t-1;
               tmp(S(uCells(c)).PixelIdxList) = 1;
               nC = nC + 1;
           end
        end
%         rgb(:,:,3) = tmp;
%         imshow(1.0*rgb)
%         pause
    end
    
end

