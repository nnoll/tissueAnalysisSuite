function [ embed ] = embeddingGrid( time_series )
    % EMBEDDING GRID 

    embed = zeros(2050,1738,3);
    [X,Y] = meshgrid(1:218,1:257);
    [Xi,Yi] = meshgrid(linspace(1,218,1738),linspace(1,257,2050));

    embed(:,:,1) = interp2(X,Y,time_series(1).eGrids{1},Xi,Yi);
    embed(:,:,2) = interp2(X,Y,time_series(1).eGrids{2},Xi,Yi);
    embed(:,:,3) = interp2(X,Y,time_series(1).eGrids{3},Xi,Yi);

end

