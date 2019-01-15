% 

for t = timePts
    tmp = cStruct(t);  
    tmp = seg.findInteriorHoles(tmp,L(:,:,t));
    tmp = seg.makeConvexArray(tmp);
    tmp.labelMat = L(:,:,t);
    tmp.embed = embed;
    [ xG{t}, yG{t} ] = measure.sparseCellGrid( tmp, 20, 20, 0 );
    disp(['t = ', str(t)])
    tic
    [ PN{t}, ERes{t}, r0{t} ] = fitDual.returnGridDual( tmp, xG{t}, yG{t} );
    toc
    [ PN{t} ] = fitDual.rescaleGridDual( PN{t} );
    [ scale{t}, goodTiles{t} ] = fitDual.propagateGridScale( PN{t}, ERes{t}, 2 );
    [ tStruct(t) ] = fitDual.storeGridMech( tmp, PN(t), scale(t), goodTiles(t) );
    [ sigma{t}, xVec{t}, yVec{t} ] = continuum.takeContinuumLimit( tStruct(t), 50, scale(t), PN(t), xG(t), yG(t), ERes(t) );
    % [ sigma{t} ] = generate.symmetrizedStress( sigma{t} );
end

% for t = timePts
%     [ tStruct(t) ] = fitDual.storeGridMech( cStruct(t), PN(t), scale(t), goodTiles(t) );
%     [ tStruct(t) ] = measure.stressTensor( tStruct(t), L(:,:,t) );
% end