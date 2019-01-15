function [ Struct, theta ] = makeStruct( p, q, Struct, extCell )

    [ ~, bulk0, ext0, ~, ~ ] = fitDual.returnGraph( Struct, extCell );
    i0 = [ bulk0, ext0 ];
    
    [ d0, ~, ~, ~, ~, involvedBonds ] = fitDual.AFN.returnBonds( Struct, i0 );
    
    x = bsxfun(@rdivide,d0*bsxfun(@times,p,q),d0*p);
    
    R = zeros(size(d0,1),1);
    Qsq = sum((d0*q).^2,2);
    Rflat = zeros(size(d0,1),1);
    
    for ii = 1:size(d0,1)
        v1 = Struct.Bdat(involvedBonds(ii)).verts(1);
        v2 = Struct.Bdat(involvedBonds(ii)).verts(2);
        r1 = double([Struct.Vdat(v1).vertxcoord,Struct.Vdat(v1).vertycoord]);
        r2 = double([Struct.Vdat(v2).vertxcoord,Struct.Vdat(v2).vertycoord]);
        Struct.Bdat(involvedBonds(ii)).rBar = x(ii,:);
        Struct.Bdat(involvedBonds(ii)).radius = mean( [sqrt(sum((r1-x(ii,:)).^2)) , sqrt(sum((r2-x(ii,:)).^2))] );
        R(ii) = Struct.Bdat(involvedBonds(ii)).radius.^2;
        Rflat(ii) = p(d0(ii,:)==1)*p(d0(ii,:)==-1)*Qsq(ii);
    end
    
    R = bsxfun(@times,R,(d0*p).^2);
    L = bsxfun(@times,d0*p,-d0);
    L = [L;ones(1,size(L,2))];
    b = [(R - Rflat);0];
    theta = L \ b;
    
end

