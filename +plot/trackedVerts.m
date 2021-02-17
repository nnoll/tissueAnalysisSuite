function [ rgb ] = trackedVerts( Struct, vertexPairs, cellPairs, skel, T )
%PLOT_TRACKED_VERTS Plots the tracer line between tracked verts
%Optional file header info (to give more details about the function than in the H1 line)
%
% Syntax:  [rgb] = trackedVerts(Struct, vertexPairs, cellPairs, skel, T)
%
% Inputs:
%    Struct - Description
%    vertexPairs - Description
%    cellPairs - Description
%    skel - 
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
% See also: TRACKEDCELLS,  TRACKEDBONDS


    rgb(:,:,1) = skel(:,:,T(1));
    rgb(:,:,3) = skel(:,:,T(2));

    tmp = zeros(size(rgb,1),size(rgb,2));
    vTrack = track.propagateVertex(vertexPairs,cellPairs,Struct,T);
    
    for v = 1:length(vTrack)
        if (vTrack(v) > 0)
            v1 = vTrack(v);
            v2 = v;
            [x, y] = plot.bresenham(Struct(T(1)).Vdat(v1).vertxcoord,Struct(T(1)).Vdat(v1).vertycoord,...
                     Struct(T(2)).Vdat(v2).vertxcoord,Struct(T(2)).Vdat(v2).vertycoord);
            ind = y + size(rgb,1)*(x-1);
            tmp(ind) = 1;
        end
    end

    tmp = imdilate(tmp,strel('disk',2));
    rgb(:,:,2) = tmp;

end

