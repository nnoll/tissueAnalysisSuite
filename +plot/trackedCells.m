function [ rgb ] = trackedCells( Struct, cellPairs, skel, T )
%PLOT_TRACKED_VERTS Plots the tracer line between tracked verts
%Optional file header info (to give more details about the function than in the H1 line)
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    Struct - Description
%    cellPairs - Description
%    skel - Description
%    T - 
%
% Outputs:
%    rgb - Description
%
% Example: 
%    Line 1 of example
%    Line 2 of example
%    Line 3 of example
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% See also: TRACKEDBONDS,  TRACKEDVERTS


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

