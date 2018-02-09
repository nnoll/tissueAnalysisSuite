function [ img2 ] = uniformBG( img, delta )
    % UNIFORM BG 

    img2 = zeros(size(img));
    S = strel('disk',30);
    G = fspecial('gaussian',40,40);
    [X,Y] = meshgrid(1:size(img,2),1:size(img,1));
    x = [X(:),Y(:)];
    
    for t = 1:size(img,3)
        t
        maxA = prctile(reshape(double(img(:,:,t)),[size(img,1)*size(img,2),1]),98);
        primeImg = mat2gray(img(:,:,t)) + delta;
        imgTmp = imfilter(primeImg,G,'replicate');
        img3 = imfilter(imopen(imgTmp,S),G,'replicate');
        
        p = load.polyfitn(x,img3(:),5);
        zG = load.polyvaln(p,x);
        zG = reshape(zG,size(img,1),size(img,2));
        img2(:,:,t) = double(max(zG(:))*primeImg./zG);
        img2(:,:,t) = maxA * img2(:,:,t) / max(max(img2(:,:,t)));

    end

    img2 = uint8(img2);
end

