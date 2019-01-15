function [ mov ] = dividingCell( L, mitCell, ID, mode )
%DIVIDINGCELL Summary of this function goes here
%   Detailed explanation goes here

    mov = plot.trackCell(mitCell(ID).label,L,mitCell(ID).t0);
    for ii = 1:size(mov,4)
        t = mitCell(ID).t0 + (ii - size(mov,4));
        S = regionprops(L(:,:,t),'Centroid');
        
        x0 = S(mitCell(ID).label(ii)).Centroid;
        if (mode == 0)
            r1 = x0 + 20*[cos(mitCell(ID).divAxis),sin(mitCell(ID).divAxis)];
            r2 = x0 - 20*[cos(mitCell(ID).divAxis),sin(mitCell(ID).divAxis)];

            mov(:,:,:,ii) = insertShape(mov(:,:,:,ii),'Line',[r1(1) r1(2) r2(1) r2(2)],'Color','blue');
            
            mov(:,:,3,ii) = mov(:,:,3,ii) - mov(:,:,2,ii);
            mov(:,:,3,ii) = imdilate(mov(:,:,3,ii),strel('disk',2));
            mov(:,:,3,ii) = mov(:,:,3,ii) + mov(:,:,2,ii);
        else
            r1 = x0 + 20*[cos(mitCell(ID).stressAxis(ii)),sin(mitCell(ID).stressAxis(ii))];
            r2 = x0 - 20*[cos(mitCell(ID).stressAxis(ii)),sin(mitCell(ID).stressAxis(ii))];

            mov(:,:,:,ii) = insertShape(mov(:,:,:,ii),'Line',[r1(1) r1(2) r2(1) r2(2)],'Color','blue');
            
            mov(:,:,3,ii) = mov(:,:,3,ii) - mov(:,:,2,ii);
            mov(:,:,3,ii) = imdilate(mov(:,:,3,ii),strel('disk',2));
            mov(:,:,3,ii) = mov(:,:,3,ii) + mov(:,:,2,ii);
        end
        
    end
    
end

