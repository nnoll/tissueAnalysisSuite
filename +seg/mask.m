function [ L ] = mask( mask )

    for t = 1:size(mask,4)
        pw = 1.0*bwmorph(mask(:,:,1,t),'diag');
        L(:,:,t) = watershed(pw);
    end

end

