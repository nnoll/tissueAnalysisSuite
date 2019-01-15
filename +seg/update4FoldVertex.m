function [ Struct ] = update4FoldVertex( Struct, PN )
    % UPDATE 4FOLD VERTEX 

    for t = 1:length(Struct)
        rPred = PN{t}.computePrimalVerts();
        [ ~, ~, ~, i0 ] = fitDual.returnGraph( Struct(t), 1 );
        for v = 1:length(Struct(t).Vdat)
            if (Struct(t).Vdat(v).old_fourfold == 1);
                vInd = (i0==v);
                Struct(t).Vdat(v).vertxcoord = rPred(vInd,1);
                Struct(t).Vdat(v).vertycoord = rPred(vInd,2);
            end
        end
    end

end

