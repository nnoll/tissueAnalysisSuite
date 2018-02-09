function [ faces ] = computeFaces( this, r )
    %COMPUTE TENSION DUAL 

    d0 = this.d0;
    d1 = this.d1;
    
    maxNVerts = max(sum(abs(d0),1));
    faces = zeros(size(d0,2),maxNVerts);
    for f = 1:size(faces,1)
        faceEdges = d0(:,f)~=0;
        faceVerts = find(sum(abs(d1(:,faceEdges)),2));
        R = bsxfun(@minus,r(faceVerts,:),mean(r(faceVerts,:),1));
        [~,ind] = sort(atan2(R(:,2),R(:,1)));
        faceVerts = faceVerts(ind);
        faces(f,1:length(faceVerts)) = faceVerts;
    end
    faces(faces==0) = nan;
    
end

