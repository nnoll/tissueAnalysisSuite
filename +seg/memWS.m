function [ L ] = memWS( mem, cellSize, strelRadius, gaussKernel, heightMiminum )
    % Use watershed method to construct the membrane segmentation.
    % 
    % Parameters
    % ----------
    % mem : membranes as images
    % cellSize : int (default=200)
    %   kernel size for laplacian of Gaussian : set to scale of curvature
    %   picking out, around a cell size or higher, in units of area (pix^2)
    % strelRadius : int (default=1)
    %   strel disk radius for dilation of segmented image
    % gaussKernel : float (default=2)
    %   kernel size for Gaussian filter. Set to a couple pixels. 
    %   Has units of length (pix).
    % heighMinimum : float (default=3.5)
    %   height of any local minima to merge, to reduce noise at rugged
    %   minima
    %
    % Returns
    % -------
    % L : label matrix
    %   edges will have zero pixel value. Cells will have integer label.
    %   The ordering of cells goes down then back up (and then back down) 
    %   as we move to the right.
    %
    % Written Nick Noll 2017-2019, Annotated NPMitchell 2019
    
    if (nargin == 1)
        cellSize = 200;
        strelRadius = 1;
        gaussKernel = 2;
        heightMiminum = 3.5;
    end
    
    % Filter output: laplacian of Gaussian sharpens, then Gaussian smooths
    h1 = fspecial('log', cellSize);
    seD1 = strel('disk', strelRadius);
    g = fspecial('gaussian', gaussKernel);
    
    % Pack label image L with watershed results
    L = zeros(size(mem));
    for t = 1:size(mem,3)
        disp(['memWS: segmenting timepoint ', num2str(t)]) 
        mem(:,:,t) = imfilter(mem(:,:,t),g);
        cyto = 1 - mem(:,:,t);
        
        mem(:,:,t) = imfilter(mem(:,:,t),h1);
        % mem(:,:,t) = imclose(mem(:,:,t),strel('disk',3));
        
        lev = graythresh(cyto);
        seed = im2bw(cyto, lev);  % consider swapping to imbinarize
        seed = imdilate(seed, seD1);
        seed = bwareaopen(seed, 25);
        
        pre_water = imhmin(mem(:,:,t), heightMiminum);
        pre_water = imimposemin(pre_water, seed);
        L(:,:,t) = watershed(pre_water);
    end


end

