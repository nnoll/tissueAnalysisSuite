delta = 1:20;
% alpha = linspace(0,20,25);
alpha = 0;
% alpha = 0;
clear LAm LAi SAm SAi
clear pStressAngle stressAngle longAngle oStressAngle tStressAngle

for jj = 1:length(alpha)

    [ Struct ] = smooth.stressTensor( Struct, L, alpha(jj) );
    [ mStruct ] = smooth.stressTensor( mStruct, L, alpha(jj) );
    
    [ Struct ] = interpolate.stressOnBndry( Struct );
    [ mStruct ] = interpolate.stressOnBndry( mStruct );

    [ mitCells ] = measure.dividingCellData( mitCells, Struct, L );
    [ mitCellsOld ] = measure.dividingCellData( mitCellsOld, mStruct, L );

    for d = delta
        ii = 1;
        for c = 1:length(mitCells)
            if ( length(mitCells(c).stressAxis) >= d );
                longAngle{jj,d}(ii) = (cos(mitCells(c).longAxis(end-d+1) - mitCells(c).divAxis)).^2;
                pStressAngle{jj,d}(ii) = (cos(mitCells(c).stressAxis(end-d+1) - mitCells(c).divAxis)).^2;
                oStressAngle{jj,d}(ii) = (cos(mitCellsOld(c).stressAxis(end-d+1) - mitCellsOld(c).divAxis)).^2;
                ii = ii + 1;
            end 
        end
    end
    
%     LAm(jj) = median(longAngle);
%     LAi(jj) = iqr(longAngle);
%     SAm(jj) = median(pStressAngle);
%     SAi(jj) = iqr(pStressAngle);
%     SAm(jj) = median(oStressAngle);
%     SAi(jj) = iqr(oStressAngle);
end

for ii = 1:size(longAngle,1)
    for t = 1:size(longAngle,2)
        diffNew(ii,t) = median(longAngle{ii,t})-median(pStressAngle{ii,t});
        diffOld(ii,t) = median(longAngle{ii,t})-median(oStressAngle{ii,t});
    end
end