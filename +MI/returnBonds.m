function [ dV, dC, t1, t2, YL, flatEdges, rho, bCells, Rpos, Rneg, S, radius ] = returnBonds( Struct )
    % RETURN BONDS 
    
    [ dC, dV, ~, iVerts, iCells ] = fitDual.ATN.computeDiffOperators( Struct, 1 );  

    [ bInd ] = generate.bondMap( Struct );
    bInd = bInd{1};
    
    t1 = zeros(size(dV,1),2);
    t2 = zeros(size(dV,1),2);
    
    rv = zeros(length(iVerts),2);
    for ii = 1:length(iVerts)
        rv(ii,1) = double(Struct.Vdat(iVerts(ii)).vertxcoord);
        rv(ii,2) = double(Struct.Vdat(iVerts(ii)).vertycoord);
    end
    
    rc = zeros(length(iCells),2);
    for ii = 1:length(iCells)
        rc(ii,1) = double(Struct.Vdat(iCells(ii)).vertxcoord);
        rc(ii,2) = double(Struct.Vdat(iCells(ii)).vertycoord);
    end
    
    rB = dV * rv;
%     cB = dC * rc;
    rho = zeros(size(dV,1),2);
%     S = sign(dot(rB,cB,2));
    YL = dC;
    bCells = zeros(size(dV,1),2);
    
    Rpos = zeros(size(dV,1),2);
    Rneg = zeros(size(dV,1),2);
    S = zeros(size(dV,1),1);
    
    for b = 1:size(dV,1)

        bVerts = [find(dV(b,:)==1),find(dV(b,:)==-1)];
        r1 = rv(bVerts(1),:);
        r2 = rv(bVerts(2),:);
        
        Rpos(b,:) = r1;
        Rneg(b,:) = r2;
        
        bCells(b,1) = find(dC(b,:)==1);
        bCells(b,2) = find(dC(b,:)==-1);
        
        flatEdges = zeros(size(dC,1),1);
        
        if (bInd(b) > 0 && Struct.Bdat(bInd(b)).radius < inf) % Finite curvature.
            
            rBar = Struct.Bdat(bInd(b)).rBar;
            if (size(rBar,1) == 2)
                rBar = rBar';
            end
            rho(b,:) = rBar;
            radius(b) = Struct.Bdat(bInd(b)).radius;
            tmp1 = r1 - rBar;
            tmp2 = r2 - rBar;
            tmp1 = tmp1' / (sqrt(sum(tmp1.^2)));
            tmp2 = tmp2' / (sqrt(sum(tmp2.^2)));
            
%             d1 = sqrt(sum( (rBar - rc(dC(b,:)==1,:)).^2 ));
%             d2 = sqrt(sum( (rBar - rc(dC(b,:)==-1,:)).^2 ));

            YL(b,:) = Struct.Bdat(bInd(b)).radius .* YL(b,:);

            if (sign(tmp1(1)*tmp2(2) - tmp1(2)*tmp2(1)) > 0)
                t1(b,:) =  [0,-1;1,0] * tmp1;
                t2(b,:) =  [0,1;-1,0] * tmp2;
                S(b) = 1;
            else
                t1(b,:) =  [0,1;-1,0] * tmp1;
                t2(b,:) =  [0,-1;1,0] * tmp2;
                S(b) = -1;
            end

        else % Flat edge
            flatEdges(b) = 1;
            t1(b,:) = -rB(b,:)/sqrt(sum(rB(b,:).^2));
            t2(b,:) = rB(b,:)/sqrt(sum(rB(b,:).^2));
            S(b) = 1;
        end            
    end   
   
    
end

