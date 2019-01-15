function [ Struct ] = stressOnBndry( Struct )
    % STRESS ON BNDRY 
    
    for t = 1:length(Struct)
        
        [ ~, ~, ~, ~, iCells ] = fitDual.ATN.computeDiffOperators( Struct(t), 1 );  
        [ ~, bulkCells ] = fitDual.returnGraph( Struct(t), 1 );
        bndryCells = iCells(~ismember(iCells,bulkCells));
        
        Rc = zeros(length(bulkCells),2);
        stress = zeros(2,2,length(bulkCells));
    
        for c = 1:length(bulkCells)
            Rc(c,:) = Struct(t).Cdat(bulkCells(c)).centroid.coord;
            stress(:,:,c) = Struct(t).Cdat(bulkCells(c)).smooth_stress;
        end
        
        RI = zeros(length(bndryCells),2);
        for c = 1:length(bndryCells)
            RI(c,:) = Struct(t).Cdat(bndryCells(c)).centroid.coord;
        end
        
        stressI = zeros(2,2,length(bndryCells));
        for ii = 1:2
            for jj = 1:2
                F = scatteredInterpolant(Rc(:,1),Rc(:,2),squeeze(stress(ii,jj,:)),'natural','linear');
                stressI(ii,jj,:) = F(RI(:,1),RI(:,2));
            end
        end

        for c = 1:length(bndryCells)
            Struct(t).Cdat(bndryCells(c)).stress = stressI(:,:,c);
            Struct(t).Cdat(bndryCells(c)).smooth_stress = stressI(:,:,c);
        end
        
    end
    
end

