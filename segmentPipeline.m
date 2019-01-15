
% membrane watershed > memWS
% L is the segmentation
%%
L = seg.memWS(mem, 200, 0, 1, 3.5);
disp('done')
%% View first and second time points to check
fig = imshow(L(:,:,1));
imshow(L(:,:,2));
disp('done')
%%
[L, Struct] = seg.generate_structs(L, 0, 1);
L = seg.removeBadCells(Struct, L);
% Once we've removed bad cells, relabel them
L = seg.relabelL(L);
[L, Struct] = seg.generate_structs(L,0,0);
disp('done')
%%
tic
% Which cells have 3 fold vertices
Struct = seg.threefold_cell(Struct);
Struct = seg.recordBonds(Struct, L);
disp('done')
toc
%%
% Compute curvature of each edge > gives centroid of circle and radii
tic
Struct = seg.curvature(Struct, size(L));
toc
%%
% T1 transition on the 4fold vertices
tic
Struct = seg.removeFourFold(Struct, size(L));
disp('done')
toc
%%
Struct = seg.makeConvexArray(Struct);
%%
tic
% Second argument is the relative coefficient of the movement of vertices
cStruct = isogonal.imposeComptCond(Struct, .5);
disp('done')
toc
%%
cStruct = seg.makeConvexArray(cStruct);


