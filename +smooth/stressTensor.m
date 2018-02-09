function [ Struct ] = stressTensor( Struct, L, smoothSize )

    for t = 1:length(Struct)
%         [ ~, ~, ~, ~, iCells ] = fitDual.ATN.computeDiffOperators( Struct(t), 1 );  
        [ ~, iCells ] = fitDual.returnGraph( Struct(t), 1 );
        for c = 1:length(iCells)
            Rc(c,:) = Struct(t).Cdat(iCells(c)).centroid.coord;
            stress(:,:,c) = Struct(t).Cdat(iCells(c)).stress;
        end
%         sum(isinf(stress(:)))
%         stress(isnan(stress)) = 0;
%         stress(isinf(stress)) = 0;
        
        S = regionprops(L(:,:,t),'Area');
        A = [S(2:end).Area];
        R0 = .5*mean(sqrt(A/pi));
        if (nargin > 2 && smoothSize > 0)
            Dc = pdist2(Rc,Rc);
            SmK = exp(-Dc.^2/(2*(smoothSize*R0)^2));
            SmK = bsxfun(@rdivide,SmK,sum(SmK,2));
            stress(1,1,:) = SmK*squeeze(stress(1,1,:));
            stress(1,2,:) = SmK*squeeze(stress(1,2,:));
            stress(2,1,:) = SmK*squeeze(stress(2,1,:));
            stress(2,2,:) = SmK*squeeze(stress(2,2,:));
        end

        for c = 1:length(iCells)
            Struct(t).Cdat(iCells(c)).smooth_stress = stress(:,:,c);
        end
    end
    
end

