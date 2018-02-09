function [ Struct ] = storeGridMech( Struct, PN, scale, goodTiles )
    
    % STORE GRID MECH 
    for t = 1:length(Struct)

        [ d0, ~, ~, ~, i0 ] = fitDual.ATN.computeDiffOperators( Struct(t), 1 );
        [ bCells ] = generate.bCells( d0, i0 );
        [ i1 ] = generate.bondMap( Struct(t) );

        T = zeros(size(bCells,1),1);
        nB = zeros(size(bCells,1),1);
        for ii = 1:size(PN{t},1)
            for jj = 1:size(PN{t},2)
                if (ismember( ii + size(PN{t},1)*(jj-1) , goodTiles{t}))
                    Tgrid = scale{t}(ii,jj) * PN{t}{ii,jj}.returnTension;
                    [ bGCells ] = generate.bCells( PN{t}{ii,jj}.d0, PN{t}{ii,jj}.cellLabels );
                    bInd = zeros(size(bGCells,1),1);
                    for b = 1:size(bGCells,1)
                        bInd(b) = find(ismember(bCells,bGCells(b,:),'rows'));
                    end
%                     bInd = ismember(bCells,bGCells,'rows');
%                     Ta = [Struct.Bdat(i1{1}(bInd)).actual_tension]';
 
%                     c(ii,jj) = corr(Tgrid,Ta);
                    T(bInd) = T(bInd) + Tgrid;
                    nB(bInd) = nB(bInd) + 1;
                end
            end
        end
%         imagesc(c)
%         pause
        T = T./nB;
        T(nB==0) = 0;

        for b = 1:length(i1{1})
            if (i1{1}(b) > 0)
                Struct(t).Bdat(i1{1}(b)).tension = T(b);
            end
        end

        for c = 1:length(i0)
            Struct(t).Cdat(i0(c)).pressure = 1;
        end

    end

end

