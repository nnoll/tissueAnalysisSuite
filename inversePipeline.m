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


%% Load in h5 file sequence in a master directory (from ilastik output).
mode = 0; % Toggle for illastik version control
Folder = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/data/48Ygal4UasCAAXmCherry/201902041850_slowdevelop_bleached_obis2_diedeventually/h5/';
imfn = '/Users/npmitchell/Dropbox/Soft_Matter/UCSB/gut_morphogenesis/data/48Ygal4UasCAAXmCherry/201902041850_slowdevelop_bleached_obis2_diedeventually/201902041850_slowdevelop_bleached_obis2_diedeventually_Cyl1_1_000011_c1_slice27.png' ;
[ mem ] = load.ilastikh5( Folder, mode );
raw = imread( imfn );

%% Segment the membrane.
% L is the label matrix
disp('segmenting the images...')
L = seg.memWS(mem, 50, 0, 1, 3.5);
% Set bond=0 and clear_border = 1
[L, Struct] = seg.generate_structs(L, 0, 1, 0, very_far);
disp('done with initial segmentation')
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
% generate the Bdat structure in Struct
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
% 'atn': Tension network inference
% 'ptn': Pressure(+Tension) network inference
atn_ptn = 'ptn'; 
extCell = 1; % The label for the external cell (surrounding void). 
             % This should honestly not be a parameter.

[PN, ERes, r0] = fitDual.returnDual(Struct, all(atn_ptn=='ptn') + 1, extCell);
disp('done with inverse')

%% Store mechanics in data structure
for t = 1:size(PN)
    % uploadMechanics is a method in the pressure.net class
    % It stores a pressure for each cdat and tension for each bdat
    [Struct(t), found_bonds] = PN{t}.uploadMechanics(Struct(t));
    Struct(t).ERes = ERes(t);
    Struct(t).PN = PN(t);
end
disp('Done storing mechanics')


%% Save to disk
outfn = split(Folder, 'h5/'); 
outfn = [ outfn{1} 'tension_net'] ;
save(outfn, 'Struct', 'ERes', 'PN')

%% Compute stress Tensor from PN
L = Struct.labelMat ;
mode = 0 ;
[ Struct ] = measure.stressTensor( Struct, L, mode ) ;

%% Plot the segmentation
L = Struct.labelMat ;
alpha = 4. ;
rgb = plot.segmentation( raw, L, alpha ) ;
if length(rgb(1,1,1,:)) > 1
    implay(rgb)
else
    imshow(rgb)
end

%% Plot the tension
imshow(L)
hold on;
% below, mode : (0 or 1) If zero, plots Struct.Bdat.tension, but if nonzero plots Struct.Bdat.actual_tension
if all(atn_ptn == 'atn')
    plot.tension(Struct, 0)
else
    plot.curvedTension(Struct, 0)
end

%% Plot the stress tensor
smoothSize = 10 ;
for t = 1:size(PN)
    plot.stressTensor( Struct(t), Struct(t).labelMat, smoothSize )
end
