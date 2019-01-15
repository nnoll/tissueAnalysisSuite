function [ Ov ] = tensionMyoTensorOvlap( Struct, L )
%TENSIONMYOTENSOROVLAP Summary of this function goes here
%   Detailed explanation goes here

    for t = 1:length(Struct)
        [ ~, ~, ~, ~, iCells ] = fitDual.ATN.computeDiffOperators( Struct(t), 1 );  

        stress = zeros(2,2,length(iCells));
        myo = zeros(2,2,length(iCells));
        
        S = regionprops(L(:,:,t),'Area');
        A = [S(iCells).Area];
        
        for c = 1:length(iCells)
            myo(:,:,c) = Struct(t).Cdat(iCells(c)).myo;
            stress(:,:,c) = Struct(t).Cdat(iCells(c)).stress;
        end
        
        lambda = mean(sum(sum(myo.*stress,2),1))/mean(sum(sum(stress.*stress,2),1));
        stress = lambda*stress;
        
        diff = myo-stress;
        delta = squeeze(sum(sum(diff.*diff,2),1)) / ...
                (mean( sqrt(sum(sum(myo.*myo,2),1)))  .* mean( sqrt(sum(sum(stress.*stress,2),1))) );
            
        Ov(t) = mean(delta);
    end
end

