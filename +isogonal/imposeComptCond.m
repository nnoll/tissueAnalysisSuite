function [ Struct, bulkVerts, bulkCells ] = imposeComptCond( Struct, kappa )
    %IMPOSE COMPT COND 

    for t = 1:length(Struct)
        t
        [ tri, bulkCells, extCells, bulkVerts, ~, rV ] = fitDual.returnGraph( Struct(t), 1 ); 
        [ extVerts, rVext ] = isogonal.returnExtVerts( bulkVerts, Struct(t) );
        [ compt, B1, B2 ] = isogonal.constructComptMatrices( tri, bulkCells, extCells, bulkVerts, extVerts, Struct(t) );
      
        x0 = double(rV(:));

        optimset = optimoptions('fminunc','Display','none','MaxFunEvals',1e6,'MaxIter',5e2,'Algorithm','quasi-newton','GradObj','on'); %,'FinDiffType','central','DerivativeCheck','on');
        x = fminunc(@(x) isogonal.comptEnergy(x, rV, rVext, kappa, sparse(compt), sparse(B1), sparse(B2) ),x0,optimset);
     
        r = reshape(x,length(x)/2,2);
        for ii = 1:length(bulkVerts)
            Struct(t).Vdat(bulkVerts(ii)).vertxcoord = double(r(ii,1));
            Struct(t).Vdat(bulkVerts(ii)).vertycoord = double(r(ii,2));
        end
    end
    
end

