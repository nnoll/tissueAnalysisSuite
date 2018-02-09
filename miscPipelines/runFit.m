N = 11;
alpha = .03;
scale = 200;
mode = 3;
gamma = .1;
beta = 1;


contin = 1;
while (contin)
    [ q, ~, b0 ] = generate.meshGrid(N,alpha);
    q(b0,:) = 1.05*q(b0,:);
    p = ones(size(q,1),1);
    z = beta*( rand(size(p)) - .5 ); 
    p = p + z;
    theta = 0 * ( rand(size(p))) ;

    [ tri, L, q, p, theta, contin ] = pressure.calculateTri( q, p, theta, b0, scale );
    clear z b0
end

PN = pressure.net(q,theta,p,tri);
T = PN.returnTension();

Struct = PN.returnStruct(gamma);
Struct = seg.threefold_cell(Struct);
%     cStruct = isogonal.imposeComptCond(Struct,.1);

[ q, p, tri, bulk0, ext0, ERes ] = fitDual.AFN.returnReducedDual( Struct, 0 );
max(sum( (PN.q([bulk0,ext0],:) - q).^2,2))
[ qI, thetaI, pI, triI, i0, ERes ] = fitDual.returnDual( Struct, mode, 0 );

if (all(~isnan(qI(:))))

    PN_i = pressure.net(qI,thetaI,pI,triI);

    [ ~, ~, ~, ~, ~, rV ] = fitDual.returnGraph( Struct, 0 );

    Ti = PN_i.returnTension();

    if (mode > 1)
        Mi = [Ti;pI/mean(pI)];
    else
        Mi = Ti;
    end

    [ bInd ] = pressure.matchBonds( PN, i0, PN_i.d0  );

    if (mode > 1)
        M = [T(bInd);p(i0)/mean(p(i0))];
    else
        M = T(bInd);
    end

    [Tm,pM] = MI.invertMech(Struct,0,mode);

    if (mode > 1)
        Mm = [Tm;pM/mean(pM)];
    else
        Mm = Tm;
    end

    if (length(Mm) > length(M))
        cOld = corr(Mm(1:length(M)),M);
    else
        cOld = corr(Mm,M(1:length(Mm)));
    end
    cNew = corr(Mi,M);
else
    cNew = 0;
    cOld = 0;
end
