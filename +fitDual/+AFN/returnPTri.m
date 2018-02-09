function [ pTri, tri ] = returnPTri( q, p, Tri )
    %RETURN PTRI 
    
    if (nargout == 1)
        
        pQ = bsxfun(@times,p,q);
        pTri = zeros(size(Tri));    
        deltaR = cat(3,pQ(Tri(:,1),:),pQ(Tri(:,2),:),pQ(Tri(:,3),:));
        deltaR = permute(bsxfun(@minus,deltaR,mean(deltaR,3)),[1,3,2]);

        angle = atan2(deltaR(:,:,2),deltaR(:,:,1));
        angle = mod(angle,2*pi);
        [~,ind] = sort(angle,2);

        pTri(:,1) = Tri((1:size(Tri,1)) + size(Tri,1)*(ind(:,1)'-1));
        pTri(:,2) = Tri((1:size(Tri,1)) + size(Tri,1)*(ind(:,2)'-1));
        pTri(:,3) = Tri((1:size(Tri,1)) + size(Tri,1)*(ind(:,3)'-1));
        
    else
        
        pQ = bsxfun(@times,p,q);
        pTri = zeros(size(Tri));    
        deltaR = cat(3,pQ(Tri(:,1),:),pQ(Tri(:,2),:),pQ(Tri(:,3),:));
        deltaR = permute(bsxfun(@minus,deltaR,mean(deltaR,3)),[1,3,2]);

        angle = atan2(deltaR(:,:,2),deltaR(:,:,1));
        angle = mod(angle,2*pi);
        [~,ind] = sort(angle,2);

        pTri(:,1) = Tri((1:size(Tri,1)) + size(Tri,1)*(ind(:,1)'-1));
        pTri(:,2) = Tri((1:size(Tri,1)) + size(Tri,1)*(ind(:,2)'-1));
        pTri(:,3) = Tri((1:size(Tri,1)) + size(Tri,1)*(ind(:,3)'-1));
        [~,ind] = min(pTri,[],2);
        pTri(ind==2,:) = circshift(pTri(ind==2,:),[0,-1]);
        pTri(ind==3,:) = circshift(pTri(ind==3,:),[0,-2]);

        tri = zeros(size(Tri));    
        deltaR = cat(3,q(Tri(:,1),:),q(Tri(:,2),:),q(Tri(:,3),:));
        deltaR = permute(bsxfun(@minus,deltaR,mean(deltaR,3)),[1,3,2]);

        angle = atan2(deltaR(:,:,2),deltaR(:,:,1));
        angle = mod(angle,2*pi);
        [~,ind] = sort(angle,2);

        tri(:,1) = Tri((1:size(Tri,1)) + size(Tri,1)*(ind(:,1)'-1));
        tri(:,2) = Tri((1:size(Tri,1)) + size(Tri,1)*(ind(:,2)'-1));
        tri(:,3) = Tri((1:size(Tri,1)) + size(Tri,1)*(ind(:,3)'-1));
        
        [~,ind] = min(tri,[],2);
        tri(ind==2,:) = circshift(tri(ind==2,:),[0,-1]);
        tri(ind==3,:) = circshift(tri(ind==3,:),[0,-2]);
    end
    
end

