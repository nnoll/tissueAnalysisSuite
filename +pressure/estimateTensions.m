function [ T, K, i1 ] = estimateTensions( p, q, Struct )

    [ ~, bulk0, ext0, ~, ~ ] = fitDual.returnGraph( Struct, 1 );
    i0 = [ bulk0, ext0 ];
    
    [ d0, ~, ~, ~, ~, i1 ] = fitDual.AFN.returnBonds( Struct, i0 );
    
    x = d0*bsxfun(@times,p,q);
        
    T = zeros(size(d0,1),1);
    deltaP = d0*p;
    
    for ii = 1:size(d0,1)
        v1 = Struct.Bdat(i1(ii)).verts(1);
        v2 = Struct.Bdat(i1(ii)).verts(2);
        r1 = double([Struct.Vdat(v1).vertxcoord,Struct.Vdat(v1).vertycoord]);
        r2 = double([Struct.Vdat(v2).vertxcoord,Struct.Vdat(v2).vertycoord]);
        T(ii) = mean( [sqrt(sum((deltaP(ii)*r1-x(ii,:)).^2)) , sqrt(sum((deltaP(ii)*r2-x(ii,:)).^2))] );
    end
    
    K = T./abs(deltaP);
    
end

