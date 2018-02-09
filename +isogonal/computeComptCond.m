function [ chi ] = computeComptCond( Struct, mode )
    %COMPUTE COMPT 
    
    for t = 1:length(Struct)
        
        [ tri, bulkCells, extCells, bulkVerts, ~, rV ] = fitDual.returnGraph( Struct(t), mode ); 

        [ extVerts, rVext ] = isogonal.returnExtVerts( bulkVerts, Struct(t) );
        [ compt, B1, B2 ] = isogonal.constructComptMatrices( tri, bulkCells, extCells, bulkVerts, extVerts, Struct(t) );
        
        rV0 = [rV;rVext];
        
        rB1 = B1*rV0;
        rB1 = bsxfun(@rdivide,rB1,sqrt(sum(rB1.^2,2)));
        rB2 = B2*rV0;
        rB2 = bsxfun(@rdivide,rB2,sqrt(sum(rB2.^2,2)));

        chi{t} = compt*log( abs(rB1(:,1).*rB2(:,2)-rB1(:,2).*rB2(:,1)) );
%         chi = ( abs(rB1(:,1).*rB2(:,2)-rB1(:,2).*rB2(:,1)) );
    end
    
end

