% nBoxes = [10,12,14,16,18,20,22,24];
nBoxes = 10;

ii = 1;
for n = nBoxes    
    [ xG, yG ] = measure.sparseCellGrid( Struct, n, n, 0 );
    [ PN, ERes, r0 ] = fitDual.returnGridDual( Struct, xG, yG );

    [ PN ] = fitDual.rescaleGridDual( PN );
    [ scale, goodTiles ] = fitDual.propagateGridScale( PN, ERes, 1 );
    [ tStruct ] = fitDual.storeGridMech( Struct, {PN}, {scale}, {goodTiles} );

    T = [];
    Ta = [];
    for b = 1:length(Struct.Bdat)
        if (~isempty(tStruct.Bdat(b).tension))
            T = [T,tStruct.Bdat(b).tension];
            Ta = [Ta,tStruct.Bdat(b).actual_tension];
        end
    end
    c(ii) = corr(T(T~=0)',Ta(T~=0)');
    ii = ii + 1;
end