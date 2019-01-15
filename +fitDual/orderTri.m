function [ nTri ] = orderTri( q, tri )
    % ORDER TRI 
    
    deltaR = cat(3,q(tri(:,1),:),q(tri(:,2),:),q(tri(:,3),:));
    deltaR = permute(bsxfun(@minus,deltaR,mean(deltaR,3)),[1,3,2]);
    
    angle = atan2(deltaR(:,:,2),deltaR(:,:,1));
    angle = mod(angle,2*pi);
    [~,ind] = sort(angle,2);
    
    nTri(:,1) = tri((1:size(tri,1)) + size(tri,1)*(ind(:,1)'-1));
    nTri(:,2) = tri((1:size(tri,1)) + size(tri,1)*(ind(:,2)'-1));
    nTri(:,3) = tri((1:size(tri,1)) + size(tri,1)*(ind(:,3)'-1));
    
end

