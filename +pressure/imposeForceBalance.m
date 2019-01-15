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

        nB = zeros(length(involvedBonds),2);
        x = zeros(length(involvedBonds),2);
        L0 = zeros(length(involvedBonds),1);
        for ii = 1:length(involvedBonds)
            x0(ii,:) = double(Struct(t).Bdat(involvedBonds(ii)).rBar);
            v1 = Struct(t).Bdat(involvedBonds(ii)).verts(1);
            v2 = Struct(t).Bdat(involvedBonds(ii)).verts(2);
            if (length(Struct(t).Bdat(involvedBonds(ii)).pix) > z)
                z = length(Struct(t).Bdat(involvedBonds(ii)).pix);
            end
            rV1(ii,:) = double([Struct(t).Vdat(v1).vertxcoord,Struct(t).Vdat(v1).vertycoord]);
            rV2(ii,:) = double([Struct(t).Vdat(v2).vertxcoord,Struct(t).Vdat(v2).vertycoord]);
            
            nB(ii,:) = (rV1(ii,:) - rV2(ii,:))*[0,-1;1,0];
            L0(ii) = sqrt(sum(nB(ii,:).^2))/2;
            nB(ii,:) = nB(ii,:)/(2*L0(ii));
            x(ii,:) = (rV1(ii,:) + rV2(ii,:))/2;
        end

        rB = zeros(length(involvedBonds),2,z);
        y0 = zeros(length(involvedBonds),1);
        for ii = 1:length(involvedBonds)
            ind = Struct(t).Bdat(involvedBonds(ii)).pix;
            [yB,xB] = ind2sub(dim,ind);
%             xB = [double(xB);rV1(ii,1);rV2(ii,1)];
%             yB = [double(yB);rV1(ii,2);rV2(ii,2)];

            rB(ii,1,1:length(xB)) = double(xB);
            rB(ii,2,1:length(xB)) = double(yB);
            y0(ii) = seg.fitToCircle(permute(rB(ii,:,1:length(xB)),[3,2,1]),nB(ii,:)',x0(ii,:)',2*L0(ii));
        end

        optimset = optimoptions('fminunc','Display','iter','TolFun',1e-4,'MaxFunEvals',1e6,'MaxIter',1e3,'Algorithm','quasi-newton','GradObj','off');
        yF = fminunc(@(y) pressure.forceBalanceEnergy(y, nB, x, L0, conBonds, rB, kappa),y0,optimset);

        for ii = 1:length(involvedBonds)
            Struct(t).Bdat(involvedBonds(ii)).rBar = (x(ii,:) + yF(ii)*nB(ii,:))';
            Struct(t).Bdat(involvedBonds(ii)).radius = sqrt(yF(ii).^2 + L0(ii).^2);
        end
    end
    
end

