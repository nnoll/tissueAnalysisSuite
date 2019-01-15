function [ X ] = computeShift( img1, img2 )
    % COMPUTE SHIFT

    w1 = hamming(size(img1,1));
    w2 = hamming(size(img1,2));
    w = w1*w2';

    F1 = fft2(w.*img1);
    F2 = fft2(w.*img2);

    R = (F1.*conj(F2))./abs(F1.*conj(F2));
    r = continuum.imgaussian(ifft2(R),2);

    [~,ind] = max(abs(r(:)));
    
    % As we used a circular shift, we need to make sure our deltaX, deltaY
    % are wrapped properly.
    
    [X(1),X(2)] = ind2sub(size(img1),ind);
    X = X - 1;
    if (abs(X(1)-1) < abs(size(img1,1)-X(1)+1))
       X(1) = -X(1)+1;
    else
       X(1) =  size(img1,1)-X(1)+1;
    end

    if abs(X(2)-1)<abs(size(img1,2)-X(2)+1)
        X(2) = -X(2)+1;
    else
        X(2) = size(img1,2)-X(2)+1;
    end

    X = X - 1;
    X = fliplr(X);

end

