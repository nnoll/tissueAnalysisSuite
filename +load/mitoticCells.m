function [ divCells ] = mitoticCells( Folder )
    % MITOTIC CELLS 
    
    files = dir(Folder);
    divCells = [];
    for n = 1:length(files)
       if (isempty(strfind(files(n).name,'.')))
           divCells = cat(4,divCells,load.loadTiff([Folder,files(n).name,'/','dividing_cells.tif']));
       end
    end
    
    set(gcf,'units','normalized','position',[0 0 1 1])

    for t = 1:size(divCells,4)
        mask = bwlabel(divCells(:,:,1,t) == 0,4);
        S = regionprops(mask,'PixelIdxList');
        imshow(divCells(:,:,:,t));
        [x,y] = ginput();
        x = round(x); y = round(y);
        cellInd = mask(y+size(mask,1)*(x-1));
        newDivCells = divCells(:,:,3,t);
        for c = 1:length(cellInd)
            newDivCells(S(cellInd(c)).PixelIdxList) = 1 - newDivCells(S(cellInd(c)).PixelIdxList);
        end
        divCells(:,:,3,t) = newDivCells;
    end
    
end

