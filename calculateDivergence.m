for t = timePts
    i1 = generate.bondMap(tStruct(t));
    T = zeros(length(i1{1}),1);
    M = zeros(length(i1{1}),1);
    for b = 1:length(i1{1})
        if (i1{1}(b) > 0)
            T(b) = tStruct(t).Bdat(i1{1}(b)).tension;
            M(b) = tStruct(t).Bdat(i1{1}(b)).myo;
        else
            T(b) = 1;
            M(b) = 1;
        end
    end

    [dC,dV,~,i0] = fitDual.ATN.computeDiffOperators(tStruct(t),1);
    rV = [[tStruct(t).Vdat(i0).vertxcoord];[tStruct(t).Vdat(i0).vertycoord]]';
    rB = dV*rV;
    rB = bsxfun(@rdivide,rB,sqrt(sum(rB.^2,2)));
    Tvec = bsxfun(@times,T,rB);
    Mvec = bsxfun(@times,M,rB);
    Mvec(isnan(Mvec)) = 1;
    
    Fnet = dV'*Tvec;
    Mnet = dV'*Mvec;
    
    FMag = sqrt(sum(Fnet.^2,2));
    MMag = sqrt(sum(Mnet.^2,2));

    FMag(FMag>.9) = 0;
    FMag(isnan(FMag)) = 0;
    
    % Keep only bulk verts
    [ ~, ~, ~, bulkVerts ] = fitDual.returnGraph( tStruct(t), 1 );
    FMag = FMag(1:length(bulkVerts));
    MMag = MMag(1:length(bulkVerts));
    rV = rV(1:length(bulkVerts),:);
    
    Sdiv(t) = nanmedian(FMag);
    Mdiv(t) = nanmean(MMag)/nanmedian(M);
    
    [ SdivMag{t} ] = interpolate.vertsOntoLab( 150, rV, FMag );
    [ MdivMag{t} ] = interpolate.vertsOntoLab( 150, rV, MMag );

end