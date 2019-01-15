function [ rgb ] = segmentation( raw, L, alpha )
%SEGMENTATION

    rgb = zeros(size(L,1),size(L,2),3,size(L,3));
    
    for t = 1:size(L,3)
        rgb(:,:,1,t) = alpha*mat2gray(raw(:,:,t));
        rgb(:,:,2,t) = .7*imdilate(L(:,:,t)==0,strel('disk',1));
    end
    

end

