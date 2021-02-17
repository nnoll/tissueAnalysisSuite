function [ L ] = mem( mem, lev )
% mem(mem, lev): membrane extraction from 2d time sequence of images
% mem : stack of membrane images
% lev : (float or int) the threshold level for binarizing the image. 
%
    L = zeros(size(mem));
    g = fspecial('gaussian',1);
    for t = 1:size(mem,3)
        if (lev == 0)
            thresh = graythresh(mem(:,:,t));
            bw = mem(:,:,t) >= thresh;
        else
            bw = mem(:,:,t) >= lev;
        end
        rgb(:,:,1) = 1*bw;
        rgb(:,:,2) = imfilter(bw,g);
        rgb(:,:,3) = watershed(imhmin(imfilter(1*bw,g),.5)) == 0;
        imshow(rgb)
        pause
        
        L(:,:,t) = bwlabel(bw);
    end


end

