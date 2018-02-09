function [ Struct ] = imposeForceBalance( Struct, dim, kappa )
    % IMPOSE FORCE BALANCE 
    
    for t = 1:length(Struct)
        t
        threefoldVerts = zeros(length(Struct(t).Vdat),1);
        involvedBonds = [];
        for v = 1:length(Struct(t).Vdat)
            if (length(Struct(t).Vdat(v).nverts) == 3 && all(Struct(t).Vdat(v).bond~=0)) % If it's three-fold
                % If all bonds have finite curvature
                if (all(~isinf([Struct(t).Bdat(Struct(t).Vdat(v).bond).radius])))
                    threefoldVerts(v) = 1; 
                    involvedBonds = [involvedBonds,Struct(t).Vdat(v).bond];
                end
            end
        end
        involvedBonds = unique(involvedBonds);
        threefoldVerts = find(threefoldVerts);

        conBonds = zeros(length(threefoldVerts),3);
        for ii = 1:length(threefoldVerts)
            conBonds(ii,1) = find(involvedBonds == Struct(t).Vdat(threefoldVerts(ii)).bond(1));
            conBonds(ii,2) = find(involvedBonds == Struct(t).Vdat(threefoldVerts(ii)).bond(2));
            conBonds(ii,3) = find(involvedBonds == Struct(t).Vdat(threefoldVerts(ii)).bond(3));
        end

        x0 = zeros(length(involvedBonds),2);
        z = 0;
        rV1 = zeros(length(involvedBonds),2);
        rV2 = zeros(length(involvedBonds),2);

        for ii = 1:length(involvedBonds)
            x0(ii,:) = double(Struct(t).Bdat(involvedBonds(ii)).rBar);
            v1 = Struct(t).Bdat(involvedBonds(ii)).verts(1);
            v2 = Struct(t).Bdat(involvedBonds(ii)).verts(2);
            if (length(Struct(t).Bdat(involvedBonds(ii)).pix) > z)
                z = length(Struct(t).Bdat(involvedBonds(ii)).pix);
            end
            rV1(ii,:) = double([Struct(t).Vdat(v1).vertxcoord,Struct(t).Vdat(v1).vertycoord]);
            rV2(ii,:) = double([Struct(t).Vdat(v2).vertxcoord,Struct(t).Vdat(v2).vertycoord]);
        end

        rB = zeros(length(involvedBonds),2,z+2);
        for ii = 1:length(involvedBonds)
            ind = Struct(t).Bdat(involvedBonds(ii)).pix;
            [y,x] = ind2sub(dim,ind);
            x = [double(x);rV1(ii,1);rV2(ii,1)];
            y = [double(y);rV1(ii,2);rV2(ii,2)];

            rB(ii,1,1:length(x)) = double(x);
            rB(ii,2,1:length(x)) = double(y);
        end
  
        optimset = optimoptions('fminunc','Display','none','TolFun',1e-4,'MaxFunEvals',1e6,'MaxIter',1e3,'Algorithm','quasi-newton','GradObj','on');
        x = fminunc(@(x) pressure.forceBalanceEnergy(x, conBonds, rV1, rV2, rB, kappa),x0,optimset);

        for ii = 1:length(involvedBonds)
            v1 = Struct(t).Bdat(involvedBonds(ii)).verts(1);
            v2 = Struct(t).Bdat(involvedBonds(ii)).verts(2);
            r1 = double([Struct(t).Vdat(v1).vertxcoord,Struct(t).Vdat(v1).vertycoord]);
            r2 = double([Struct(t).Vdat(v2).vertxcoord,Struct(t).Vdat(v2).vertycoord]);
            Struct(t).Bdat(involvedBonds(ii)).rBar = x(ii,:);
            Struct(t).Bdat(involvedBonds(ii)).radius = mean( [sqrt(sum((r1-x(ii,:)).^2)) , sqrt(sum((r2-x(ii,:)).^2))] );
        end
    end
    
end

