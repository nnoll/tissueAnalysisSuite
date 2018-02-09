
nBoxes = 5:20;
prctile = linspace(0,.25,3);

clear c
for jj = 1:length(prctile)
    for ii = 1:length(nBoxes)
        clear PN ERes r0

        ii
        [ xG, yG ] = measure.sparseCellGrid( mStruct, nBoxes(ii), nBoxes(ii), prctile(jj) );
        [ PN{1}, ERes{1}, r0{1} ] = fitDual.returnGridDual( mStruct, xG, yG );
        [ PN{1} ] = fitDual.rescaleGridDual( PN{1} );
        [ scale{1}, goodTiles{1} ] = fitDual.propagateGridScale2( PN{1}, ERes{1}, 1 );
        [ c(jj,ii) ] = fitDual.storeGridMech( mStruct, PN, scale, goodTiles );

    end
end