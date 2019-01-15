function [ L2 ] = relabelL( L )
    
    L2 = zeros(size(L));
    for t = 1:size(L,3)
        tmp = zeros(size(L,1),size(L,2));
        Lt = L(:,:,t);
        ids = sort(unique(Lt(:)));
        ids(ids==0) = [];
        for ii = 1:length(ids)
           tmp(Lt == ids(ii)) = ii; 
        end
        L2(:,:,t) = tmp;
    end


end

