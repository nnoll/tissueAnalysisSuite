function [ E, dE ] = energy( x, nCon, nB, D1, D2, D3, L )
%ENERGY 

    r0 = reshape(x,length(x)/2,2);
    
    v1 = (r0(nCon(:,1),:) - r0(1:nB,:));
    v2 = (r0(nCon(:,2),:) - r0(1:nB,:));
    v3 = (r0(nCon(:,3),:) - r0(1:nB,:));
    
    V1 = dot(v1,v1,2);
    V2 = dot(v2,v2,2);
    V3 = dot(v3,v3,2);
    
%     v1 = bsxfun(@rdivide,v1,sqrt(V1));
%     v2 = bsxfun(@rdivide,v2,sqrt(V2));
%     v3 = bsxfun(@rdivide,v3,sqrt(V3));

    E1 = dot(v1,v2,2);
    E2 = dot(v2,v3,2);
    E3 = dot(v3,v1,2);
    
    E = sum( (E1-D1).^2 ) + sum( (E2-D2).^2 )  + sum( (E3-D3).^2 );
    E = E / (6*size(x,1));
    
    if (nargout > 1)
       dE = zeros(size(r0));

       % Add contributions from self terms
       dE(1:nB,:) = dE(1:nB,:) + bsxfun(@times,(E1-D1),(2*r0(1:nB,:) -r0(nCon(:,1),:) - r0(nCon(:,2),:)));
       dE(1:nB,:) = dE(1:nB,:) + bsxfun(@times,(E2-D2),(2*r0(1:nB,:) -r0(nCon(:,2),:) - r0(nCon(:,3),:)));
       dE(1:nB,:) = dE(1:nB,:) + bsxfun(@times,(E3-D3),(2*r0(1:nB,:) -r0(nCon(:,3),:) - r0(nCon(:,1),:)));
       
       % Add contributions from cross terms
       edge1 = bsxfun(@times,(E1-D1), r0(nCon(:,2),:) - r0(1:nB,:)) + bsxfun(@times,(E3-D3), r0(nCon(:,3),:) - r0(1:nB,:));
       edge2 = bsxfun(@times,(E1-D1), r0(nCon(:,1),:) - r0(1:nB,:)) + bsxfun(@times,(E2-D2), r0(nCon(:,3),:) - r0(1:nB,:));
       edge3 = bsxfun(@times,(E2-D2), r0(nCon(:,2),:) - r0(1:nB,:)) + bsxfun(@times,(E3-D3), r0(nCon(:,1),:) - r0(1:nB,:));
       
       dE(:,1) = dE(:,1) + L * [edge1(:,1);edge2(:,1);edge3(:,1)];
       dE(:,2) = dE(:,2) + L * [edge1(:,2);edge2(:,2);edge3(:,2)];

       dE = dE(:)/(3*size(x,1));
    end
end

