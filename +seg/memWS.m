function [ L ] = memWS( mem, H, S, G, HM )

    if (nargin == 1)
        H = 200;
        S = 1;
        G = 2;
        HM = 3.5;
    end
    
    h1 = fspecial('log',H);
    seD1 = strel('disk',S);
    g = fspecial('gaussian',G);
    
    L = zeros(size(mem));
    for t = 1:size(mem,3)
        mem(:,:,t) = imfilter(mem(:,:,t),g);
        cyto = 1 - mem(:,:,t);
        
        mem(:,:,t) = imfilter(mem(:,:,t),h1);
%         mem(:,:,t) = imclose(mem(:,:,t),strel('disk',3));
        
        lev = graythresh(cyto);
        seed = im2bw(cyto,lev);
        seed = imdilate(seed,seD1);
        seed = bwareaopen(seed,25);
        
        pre_water = imhmin(mem(:,:,t),HM);
        pre_water = imimposemin(pre_water,seed);
        L(:,:,t) = watershed(pre_water);
    end


end

