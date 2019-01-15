function [ rgb ] = trackedCells( Struct, cellPairs, skel, T )
%PLOT_TRACKED_VERTS Plots the tracer line between tracked verts

    rgb(:,:,1) = skel(:,:,T(1));
    rgb(:,:,3) = skel(:,:,T(2));

    tmp = zeros(size(rgb,1),size(rgb,2));
    cTrack = track.propagateCells(cellPairs,T);
    
    for c = 1:length(cTrack)
        if (cTrack(c) > 0)
            c1 = c;
            c2 = cTrack(c);
            [x, y] = plot.bresenham(Struct(T(1)).Cdat(c2).centroid.coord(1),...
                     Struct(T(1)).Cdat(c2).centroid.coord(2),Struct(T(2)).Cdat(c1).centroid.coord(1),...
                     Struct(T(2)).Cdat(c1).centroid.coord(2));
            ind = y + size(rgb,1)*(x-1);
            tmp(ind) = 1;
        end
    end

    tmp = imdilate(tmp,strel('disk',2));
    rgb(:,:,2) = tmp;

end

