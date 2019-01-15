function [ cellPairs ] = cells( L, cutoff )

    S0 = regionprops(L(:,:,1),'PixelIdxList','Centroid');
    R0 = vertcat(S0.Centroid);
    
    cellPairs = cell(size(L,3),1);
    
    for t = 2:size(L,3)
       t
       S = regionprops(L(:,:,t),'PixelIdxList','Centroid'); 
       R = vertcat(S.Centroid);

       % Pre-filtering step based upon distance.
       D = pdist2(R,R0);
       [currentCells,oldCells] = find(D < cutoff);

       Ov = inf*ones(size(D));
       for ii = 1:length(oldCells)
           sharedArea = sum(ismembc(S(currentCells(ii)).PixelIdxList,S0(oldCells(ii)).PixelIdxList));

           Ov(currentCells(ii),oldCells(ii)) = ... 
           (length(S(currentCells(ii)).PixelIdxList) + length(S0(oldCells(ii)).PixelIdxList))/sharedArea - 2;
       end
       
       [ M ] = track.munkres(Ov);
       matchedCells = find(M);
       cellPairs{t} = zeros(length(S),1);
       cellPairs{t}(matchedCells) = M(matchedCells);
       
       S0 = S;
       R0 = R;
    end

end

