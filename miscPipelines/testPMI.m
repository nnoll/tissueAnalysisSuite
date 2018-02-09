N = 11;
alpha = .03;
beta = .8;
gamma = 0;
mode = 3;

[ q, ~, b0 ] = generate.meshGrid(N,alpha);
q(b0,:) = 1.05*q(b0,:);
p = ones(size(q,1),1);
z = beta*( rand(size(p)) - .5 ); 
p = p + z;
% p(b0) = 1;
theta = .2 * ( rand(size(p))) ;

[ tri, L, q, p, theta, contin ] = pressure.calculateTri( q, p, theta, b0 );
clear z b0

PN = pressure.net(q,theta,p,tri);
Struct = PN.returnStruct(gamma);
Struct = seg.threefold_cell(Struct);

[ qI, thetaI, pI, triI, i0, ERes ] = fitDual.returnDual( Struct, mode );
D = pdist2(q,qI);
ind = track.munkres(D);

if (sum(ind==0) == 0)
    qI = qI(ind,:);
    thetaI = thetaI(ind);
    pI = pI(ind);

    PN_i = pressure.net(qI,thetaI,pI,tri);
    [bInd,cInd] = MI.matchBondsAndCells(Struct,PN,0);
    [Tm,Pm] = MI.invertMech(Struct,0,mode);

    T = [PN.returnTension();p];
    Ti = [PN_i.returnTension();pI];
    Tm = [Tm(bInd);Pm(cInd)];

    deltaP = mean(abs(PN.d0*PN.p));

    bVerts = vertcat(Struct.Bdat.verts);
    clear r
    r(:,1) = [Struct.Vdat.vertxcoord];
    r(:,2) = [Struct.Vdat.vertycoord];
    rB = r(bVerts(:,1),:) - r(bVerts(:,2),:);
    rB = sqrt(sum(rB.^2,2));

    deltaR = median(rB);
    dCorr = corr(T,Ti);
    mCorr = corr(Tm,T);

    clear r rB bVerts
    
end
