function [ Struct ] = curvedVM( Struct, dim )

    % CURVED VERTEX MODEL 
    for t = 1:length(Struct)
       
        [ dC, dV, ~, iVerts, iCells ] = fitDual.ATN.computeDiffOperators( Struct(t), 1 );  
        [ involvedBonds ] = generate.bondMap(Struct(t));
        involvedBonds = involvedBonds{1};
        
        z = 0;
        bVerts = zeros(length(involvedBonds),2);
        badBonds = [];
        for b = 1:length(involvedBonds)
            if (involvedBonds(b) > 0)
                v1 = find(iVerts == Struct(t).Bdat(involvedBonds(b)).verts(1));
                v2 = find(iVerts == Struct(t).Bdat(involvedBonds(b)).verts(2));
                if (length(Struct(t).Bdat(involvedBonds(b)).pix) > z)
                    z = length(Struct(t).Bdat(involvedBonds(b)).pix);
                end
                bVerts(b,:) = [v1,v2];
            else
                badBonds = [badBonds,b];
            end
        end
        involvedBonds(badBonds) = [];
        bVerts(badBonds,:) = [];

        rV = [Struct(t).Vdat(iVerts).vertxcoord;Struct(t).Vdat(iVerts).vertycoord]';
        
        rB = zeros(length(involvedBonds),2,z+2);
        for ii = 1:length(involvedBonds)
            ind = Struct(t).Bdat(involvedBonds(ii)).pix;
            [y,x] = ind2sub(dim,ind);
            x = [double(x);rV(bVerts(ii,1),1);rV(bVerts(ii,2),1)];
            y = [double(y);rV(bVerts(ii,1),2);rV(bVerts(ii,2),2)];

            rB(ii,1,1:length(x)) = double(x);
            rB(ii,2,1:length(x)) = double(y);
        end
        
        mask = permute(rB(:,1,:) ~= 0,[1,3,2]);
        optimset = optimoptions('fminunc','Display','iter','TolFun',1e-6,'MaxFunEvals',1e6,'MaxIter',1e3,'Algorithm','quasi-newton'); %,'GradObj','on');
        x = fminunc(@(q) generate.segEnergy(q, bVerts, rB, mask),q0,optimset);
        
        for v = 1:length(iVerts)
            Struct(t).Vdat(iVerts(v)).vertxcoord = x(v,1);
            Struct(t).Vdat(iVerts(v)).vertycoord = x(v,2);
        end
        
        Struct(t) = seg.curvature(Struct(t),dim);
    end

end

