% This code will generate a data structure called Struct which contains the
% following fields:
% vdat: 
%     nverts : neighbor list for vertices
%     ncells : 
%     vertxcoord : column of data in which vertex lives
%     vertycoord : row of data in which each vertex lives
% cdat : cell data    
%     ncells : indices of neighboring cells
%     nverts : vertices that define the cell
%     centroid.coord(1) x position of cell centroid
%     centroid.coord(2) y position of cell centroid
% bdat : 
%     nverts
%     ncells
%     pix : linear indices of the pixels associated with that bond

% Note on exterior calculus objects
% ---------------------------------
% d0 and d1 are matrices that take derivatives.
% d0 is an e x c matrix of exterior derivatives with +1 and -1s
% at the endpts of each bond.
% d1 is a v x e matrix of exterior derivatives. Upstream is +1,
% downstream is -1 when moving counterclockwise around a tension plaquette.

%% Parameters
very_far = 150 ;


%% Load in h5 from ilastik.
mode = 0; % Toggle for illastik version control
Folder = 'h5/';
[ mem ] = load.ilastikh5( Folder, mode );

%% Segment the membrane.
L = seg.memWS(mem, 50, 0, 1, 3.5);
% Set bond=0 and clear_border = 1
[L, Struct] = seg.generate_structs(L, 0, 1, 0, very_far);
% Bad cells are bubble cells, which is a segmentation that forked and
% reconnected.
L = seg.removeBadCells(Struct, L);
disp('done removing bad cells')
% Now change label matrix after removing bad cells
L = seg.relabelL(L);
% Now also synchronize Struct after removing bad cells
[L,Struct] = seg.generate_structs(L, 0, 0, 0, very_far);
disp('done with segmentation')

%% Prepare data structure for inverse.
% put a parameter in the cdat of Struct, a boolean of whether every vertex
% is 3-fold.
Struct = seg.threefold_cell(Struct);
% generate the bdat structure in Struct
Struct = seg.recordBonds(Struct, L);
disp('generated the bond structure')
% Segment the curvature of each bond
Struct = seg.curvature(Struct, size(L));
disp('segmented the curvature of each bond')
% Remove all fourfold vertices, recursively if there are z>4
Struct = seg.removeFourFold(Struct);
disp('removed fourfold vertices')
% The inverse is ill-posed if we have convex cells, so hack those to be
% convex
Struct = seg.makeConvexArray(Struct);
disp('done with data preparation')

%%
% This does a soft version of compatibility constraint using Monte Carlo
% sampling to make edges lie along a line. May not be necessary.
% Struct = isogonal.imposeComptCond(Struct, .5);
% Struct = seg.makeConvexArray(Struct);

%% 
for t = 1:size(L, 3)
    Struct(t).labelMat = L(:, :, t);
end
clear L
disp('done with putting L data into Struct')

%% Invert mechanics.
mode = 3; % Pressure network inference
extCell = 1; % This should honestly not be a parameter.
[PN, ERes, r0] = fitDual.returnDual(Struct, mode, extCell);
disp('done with inverse')

%% Store mechanics in data structure
for t = 1:size(PN)
    % uploadMechanics is a method in the pressure.net class
    % It stores a pressure for each cdat and tension for each bdat
    [Struct(t), found_bonds] = PN{t}.uploadMechanics(Struct(t), ERes);
    Struct(t).ERes = ERes(t);
    Struct(t).PN = PN(t);
end
disp('Done storing mechanics')


%% Save to disk
saveas(outfn, [Struct, ERes, PN])
