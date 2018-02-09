function [ L, scale ] = weightedVoronoi( r, d, dim )
%GENERATEWEIGHTEDVORONOI 

    r = round(r);
    
    D = inf*ones(dim);
    for n = 1:size(r,2)
        bw = zeros(dim);
        bw(r(2,n),r(1,n)) = 1;
        D = min(D,bwdist(bw)/d(n));
    end
    

    L = watershed(D);
    
    %Relabel L
    labels = L(r(2,:) + size(L,1)*(r(1,:)-1));
    Ltmp = zeros(size(L));
    for ii = 1:length(labels)
       Ltmp(L == labels(ii)) = ii; 
    end
    L = Ltmp;
    
    [Tri] = delaunayTriangulation(r');
    
    % Make Adj matrix.
    Adj = zeros(size(r,2));
    for f = 1:size(Tri,1)
       Adj(Tri(f,1),Tri(f,2)) = 1; 
       Adj(Tri(f,2),Tri(f,1)) = 1; 

       Adj(Tri(f,3),Tri(f,2)) = 1; 
       Adj(Tri(f,2),Tri(f,3)) = 1; 
       
       Adj(Tri(f,1),Tri(f,3)) = 1; 
       Adj(Tri(f,3),Tri(f,1)) = 1; 
    end
    
    curvature = zeros(1,sum(Adj(:)));
    n = 1;
    for ii = 1:size(Adj,1)
       for jj = find(Adj(ii,:))
          curvature(n) = abs(1/d(ii) - 1/d(jj))./norm(r(:,ii)-r(:,jj));
          n = n + 1;
       end
    end
    scale = mean(curvature);

end

