function [ rgb ] = trackedBonds( Struct, cPair, skel, T )
%TRACKEDBONDS Plots the tracer line between tracked bonds?
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%Optional file header info (to give more details about the function than in the H1 line)
%
% Syntax:  [output1,output2] = function_name(input1,input2,input3)
%
% Inputs:
%    input1 - Description
%    input2 - Description
%    input3 - Description
%
% Outputs:
%    output1 - Description
%    output2 - Description
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
% See also: OTHER_FUNCTION_NAME1,  OTHER_FUNCTION_NAME2

    rgb(:,:,1) = skel(:,:,T(1));
    rgb(:,:,3) = skel(:,:,T(2));

    tmp = zeros(size(rgb,1),size(rgb,2));
    bTrack = track.propagateBond(Struct,cPair,T);
    
    b0Verts = [Struct(T(1)).Bdat.verts];
    bVerts = [Struct(T(2)).Bdat.verts];
    r0 = [Struct(T(1)).Vdat.vertxcoord;Struct(T(1)).Vdat.vertycoord];
    r = [Struct(T(2)).Vdat.vertxcoord;Struct(T(2)).Vdat.vertycoord];
    rB0 = .5*(r0(:,b0Verts(1,:)) + r0(:,b0Verts(2,:)));
    rB = .5*(r(:,bVerts(1,:)) + r(:,bVerts(2,:)));
    
    for b = 1:length(bTrack)
        if (bTrack(b) > 0)
            b1 = bTrack(b);
            b2 = b;
            [x, y] = plot.bresenham(rB0(1,b1),rB0(2,b1),rB(1,b2),rB(2,b2));
            ind = y + size(rgb,1)*(x-1);
            tmp(ind) = 1;
        end
    end

    tmp = imdilate(tmp,strel('disk',2));
    rgb(:,:,2) = tmp;

end

