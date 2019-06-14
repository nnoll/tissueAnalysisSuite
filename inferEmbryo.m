% Script for inferring tensions in 3D using tangent patches
% tmp.embed is a map from (u,v) in charts to (x,y,z) in R3.

for t = timePts
    tmp = cStruct(t);  
    tmp = seg.findInteriorHoles(tmp,L(:,:,t));
    tmp = seg.makeConvexArray(tmp);
    % fix the label matrix from watershed
    tmp.labelMat = L(:,:,t);
    % add attribute that is the output from imsane
    tmp.embed = embed;
    % Define a xy grid of local patches to map to tangent planes
    % Consider changing so that there are overlapping patches from a dual
    % of a Vogel disk with overlapping regions.
    % For example, tile the disk with points (Vogel), take dual to get a
    % Voronoi tesselation. Dilate each voronoi cell by 20% so they overlap.
    % Use overlap to constrain relative scale of the stresses inferred.
    [ xG{t}, yG{t} ] = measure.sparseCellGrid( tmp, 20, 20, 0 );
    disp(['t = ', str(t)])
    
    % Turn points into patches and do the mechanical inverse in each.
    tic
    [ PN{t}, ERes{t}, r0{t} ] = fitDual.returnGridDual( tmp, xG{t}, yG{t} );
    toc
    
    % Set the average tensions in each patch to 1. 
    [ PN{t} ] = fitDual.rescaleGridDual( PN{t} );
    % Rescale the tensions in each patch to match overlapping portions. 
    % goodTiles is a binary mask in the defined patches. If there a bad fit
    % where the residuals are high, goodTiles is a mask.
    [ scale{t}, goodTiles{t} ] = fitDual.propagateGridScale( PN{t}, ERes{t}, 2 );
    
    % Store the mechanics in a structure
    [ tStruct(t) ] = fitDual.storeGridMech( tmp, PN(t), scale(t), goodTiles(t) );
    % Coarse grain the stresses
    [ sigma{t}, xVec{t}, yVec{t} ] = continuum.takeContinuumLimit( tStruct(t), 50, scale(t), PN(t), xG(t), yG(t), ERes(t) );
    % [ sigma{t} ] = generate.symmetrizedStress( sigma{t} );
end

% for t = timePts
%     [ tStruct(t) ] = fitDual.storeGridMech( cStruct(t), PN(t), scale(t), goodTiles(t) );
%     [ tStruct(t) ] = measure.stressTensor( tStruct(t), L(:,:,t) );
% end