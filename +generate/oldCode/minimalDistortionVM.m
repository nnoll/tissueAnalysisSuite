function [ rv, ERes ] = minimalDistortionVM( Struct, bulkVerts, extVerts, g )
    %MINIMALDISTORTIONVM 
    
    iVerts = [bulkVerts,extVerts];
    nB = length(bulkVerts);
    
    r0 = zeros(length(iVerts),2);
    for ii = 1:length(iVerts)
        r0(ii,1) = double(Struct.Vdat(iVerts(ii)).vertxcoord);
        r0(ii,2) = double(Struct.Vdat(iVerts(ii)).vertycoord);
    end
    
    if (isa(g,'function_handle'))
        g0 = g(r0(1:nB,1),r0(1:nB,2));
    else
        r0G = round(r0);

        g0 = zeros(nB,3);
        g1 = g(:,:,1);     g2 = g(:,:,2);     g3 = g(:,:,3);
        g0(:,1) = g1(r0G(1:nB,2) + size(g,1)*(r0G(1:nB,1)-1));
        g0(:,2) = g2(r0G(1:nB,2) + size(g,1)*(r0G(1:nB,1)-1));
        g0(:,3) = g3(r0G(1:nB,2) + size(g,1)*(r0G(1:nB,1)-1));
    end
    
    nCon = zeros(length(bulkVerts),3);
    for ii = 1:length(bulkVerts)
        nCon(ii,:) = [find(iVerts==Struct.Vdat(bulkVerts(ii)).nverts(1)),...
                      find(iVerts==Struct.Vdat(bulkVerts(ii)).nverts(2)),...
                      find(iVerts==Struct.Vdat(bulkVerts(ii)).nverts(3))];
        
    end
    
    v1 = (r0(nCon(:,1),:) - r0(1:nB,:));
    v2 = (r0(nCon(:,2),:) - r0(1:nB,:));
    v3 = (r0(nCon(:,3),:) - r0(1:nB,:));

    V1 = v1(:,1) .* g0(:,1) .* v1(:,1) + v1(:,1) .* g0(:,2) .* v1(:,2) + ...
         v1(:,2) .* g0(:,2) .* v1(:,1) + v1(:,2) .* g0(:,3) .* v1(:,2);
     
    V2 = v2(:,1) .* g0(:,1) .* v2(:,1) + v2(:,1) .* g0(:,2) .* v2(:,2) + ...
         v2(:,2) .* g0(:,2) .* v2(:,1) + v2(:,2) .* g0(:,3) .* v2(:,2);
     
    V3 = v3(:,1) .* g0(:,1) .* v3(:,1) + v3(:,1) .* g0(:,2) .* v3(:,2) + ...
         v3(:,2) .* g0(:,2) .* v3(:,1) + v3(:,2) .* g0(:,3) .* v3(:,2);
     
    % Normalize the results
%     v1 = bsxfun(@rdivide,v1,sqrt(V1));
%     v2 = bsxfun(@rdivide,v2,sqrt(V2));
%     v3 = bsxfun(@rdivide,v3,sqrt(V3));
    
    % Compute inner products
    D1 = v1(:,1) .* g0(:,1) .* v2(:,1) + v1(:,1) .* g0(:,2) .* v2(:,2) + ...
         v1(:,2) .* g0(:,2) .* v2(:,1) + v1(:,2) .* g0(:,3) .* v2(:,2);
     
    D2 = v2(:,1) .* g0(:,1) .* v3(:,1) + v2(:,1) .* g0(:,2) .* v3(:,2) + ...
         v2(:,2) .* g0(:,2) .* v3(:,1) + v2(:,2) .* g0(:,3) .* v3(:,2);
     
    D3 = v3(:,1) .* g0(:,1) .* v1(:,1) + v3(:,1) .* g0(:,2) .* v1(:,2) + ...
         v3(:,2) .* g0(:,2) .* v1(:,1) + v3(:,2) .* g0(:,3) .* v1(:,2);
    
    optimset = optimoptions('fminunc','Display','none','Algorithm','quasi-newton','GradObj','on','MaxIter',1000,'MaxFunEvals',50000);
    
    L = zeros(length(iVerts),3*nB);
    L(nCon(:,1) + length(iVerts)*((1:nB)'-1)) = 1;
    L(nCon(:,2) + length(iVerts)*(((nB+1):2*nB)'-1)) = 1;
    L(nCon(:,3) + length(iVerts)*(((2*nB+1):3*nB)'-1)) = 1;

    L = sparse(L);
    energyFunc = @(x) generate.energy(x,nCon,nB,D1,D2,D3,L);

    [x,ERes] = fminunc(energyFunc,r0(:),optimset);
    rv = reshape(x,length(x)/2,2);
    
%     plot.skel(Struct,'k',0)
%     hold on
%     scatter(r0(:,1),r0(:,2),'b','filled')
%     scatter(rv(:,1),rv(:,2),'ro')
%     pause
end

