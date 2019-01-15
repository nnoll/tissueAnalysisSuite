function [ rgb ] = cellHM( Theta, bulkCells, L, Struct, T )
%PLOT_THETA Plots the fitted theta back onto the movie

%%Inputs 
%1. ThetaC - Fitted dilation
%2. three_cells - List of cell labels that were fitted
%3. mem - Membrane channel
%4. L - Watershed
%5. LUT - Cell tracker
%6. T - [Tstart, Tend]

%%Outputs
%1. rgb - output movie

    if (isempty(bulkCells))
        [ ~, ~, ~, ~, bulkCells ] = fitDual.ATN.computeDiffOperators( Struct, 1 );  
    end

    S = regionprops(L(:,:,T),'PixelIdxList');
    rgb = zeros(size(L(:,:,T)));

    n = 1;
    for c = bulkCells
        rgb(S(c).PixelIdxList) = Theta(n);
        n = n + 1;
    end

    extCells = 1:length(S);
    extCells(bulkCells) = [];

    mask = zeros(size(L(:,:,T)));
    for c = extCells
        mask(S(c).PixelIdxList) = 1;
    end
    mask = imdilate(mask,strel('disk',1));

    rgb = 1.0*rgb;
    rgb(mask==1) = 0;
    rgb(imdilate(L(:,:,T)==0,strel('disk',1))) = 0;

end

