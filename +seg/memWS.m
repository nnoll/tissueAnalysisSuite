function [ L ] = memWS( mem, H, S, G, HM )
    % Use watershed method to construct the membrane segmentation.
    % 
    % Parameters
    % ----------
    % mem : membranes
    % H : float
    %   kernel size for laplacian of Gaussian : set to scale of curvature
    %   picking out, around a cell size or higher, in units of area (pix^2)
    % G : float
    %   kernel size for Gaussian filter. Set to a couple pixels. 
    %   Has units of length (pix).
    % HM : 
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
        H = 200;
        S = 1;
        G = 2;
        HM = 3.5;
    end
    
    % Filter output: laplacian of Gaussian sharpens, then Gaussian smooths
    h1 = fspecial('log', H);
    seD1 = strel('disk', S);
    g = fspecial('gaussian', G);
    
    % Pack label image L with watershed results
    L = zeros(size(mem));
    for t = 1:size(mem,3)
        disp(['memWS: segmenting timepoint ', num2str(t)]) 
        mem(:,:,t) = imfilter(mem(:,:,t),g);
        cyto = 1 - mem(:,:,t);
        
        mem(:,:,t) = imfilter(mem(:,:,t),h1);
        % mem(:,:,t) = imclose(mem(:,:,t),strel('disk',3));
        
        lev = graythresh(cyto);
        seed = im2bw(cyto, lev);
        seed = imdilate(seed, seD1);
        seed = bwareaopen(seed, 25);
        
        pre_water = imhmin(mem(:,:,t), HM);
        pre_water = imimposemin(pre_water, seed);
        L(:,:,t) = watershed(pre_water);
    end


end

