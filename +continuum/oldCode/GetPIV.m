function [VX,VY] = GetPIV(im1,im2,X1,Y1,EdgeLength)
    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    %   GetPIV computes PIV flow field estimate based on image1 (im1) and
    %   image2 (im2) using the Phase Correlation method implemented in xcorr2fft. 
    %   im1 & im2 are assumed to have the same dimensions. The grid X1,Y1
    %   is assmued to be contained in the image domain with finite
    %   EdgeLength defining size of PIV box. 
    %   Output are components of the flow fiel on the grid
    %   
    %   Written by: Sebastian J Streichan, KITP, February 01, 2013
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    
    si      = size(im1);
    nb      = ceil(si/EdgeLength);
    counter = 0;
    shiftx  = zeros(1,nb(1)*nb(2));
    shifty  = shiftx;
    x       = shiftx;
    y       = shifty;
    for i = 1 : nb(1)
        for j = 1 : nb(2)

            counter = counter+1;
            temp1 = im1(max(1,((i-2)*EdgeLength+1)):min((i+1)*EdgeLength,si(1)),(max(1,(j-2)*EdgeLength+1)):min((j+1)*EdgeLength,si(2)));
            temp2 = im2(max(1,((i-2)*EdgeLength+1)):min((i+1)*EdgeLength,si(1)),(max(1,(j-2)*EdgeLength+1)):min((j+1)*EdgeLength,si(2)));    
            
            [shiftx(counter),shifty(counter)] = xcorr2fft(temp1,temp2);
            x(counter) = (i-.5)*EdgeLength;
            y(counter) = (j-.5)*EdgeLength;
        end
    end
    
    % Put the flow field on a grid 
    VX  = griddata(x,y,shiftx,X1,Y1);
    VY  = griddata(x,y,shifty,X1,Y1);
end