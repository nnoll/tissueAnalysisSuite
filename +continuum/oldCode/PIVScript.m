%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Using GetPIV we computes the flow field estimate based on an image sequence 
%   im1 & im2 are assumed to have the same dimensions. The grid X1,Y1
%   is assmued to be contained in the image domain with finite
%   EdgeLength defining size of PIV box. 
%   This script generates a movie displaying the original image with flow
%   field overlayed. 
%   
%   Written by: Sebastian J Streichan, KITP, February 14, 2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Name    = '../Layer05-Equalized.tif';

StackSize   = length(imfinfo(Name));
EdgeLength  = 15;   % Length of box edges in pixels; 
EdgeLength2 = 5;    % Length of boxedges in interpolated field
isf         = .4;   % image scaling factor. 
step        = 3;    % step in timeframes.
smooth      = 1;    % set to 1 if gaussian smoothing is desired
KernelSize  = 7;    % Smoothing kernel size
sigma       = 1.7;  % standard deviation of gaussian kernel
im1 = imread(Name,1);
im1 = imresize(im1,isf,'bicubic');

% define the grid on which to compute the flow field
[X1,Y1] = meshgrid(EdgeLength/2:EdgeLength:size(im1,1)-EdgeLength/2,EdgeLength/2:EdgeLength:size(im1,2)-EdgeLength/2); 

for t = 1:  StackSize-step

    % read the image and scale
    im1     = imread(Name,t);
    im2     = imread(Name,t+step);
    im1     = imresize(im1,isf,'bicubic'); % rescale image if desired
    im2     = imresize(im2,isf,'bicubic');

    % compute the piv flow field
    [VX,VY] = GetPIV(im1,im2,X1,Y1,EdgeLength); 
    
    % smooth if desired
    if smooth == 1
        VX  = imfilter(VX,fspecial('gaussian',KernelSize,sigma));
        VY  = imfilter(VY,fspecial('gaussian',KernelSize,sigma));
    end
    
    % Display image and overlay flow field.
    imshow(im1',[])
    hold on 
    quiver(X1,Y1,VX,VY,10,'g-')
    
    % records a movie
    M(t)    = getframe; 
end
%play a movie;
implay(M) 