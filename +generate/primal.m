function [ d0, d1, bndry0, bndry1, bndry2 ] = primal( Mesh )
    %GENERATE TRIANGULATION 
    
    % Return boundary verts
    bndry1 = Mesh.freeBoundary();
    bndry0 = unique(bndry1(:));
    
    % Define exterior derivatives on vertices and edges.
    edges = Mesh.edges;
    d0 = zeros(size(edges,1),size(Mesh.Points,1));
    d0((1:size(edges,1)) + size(edges,1)*(edges(:,1)'-1)) = 1;
    d0((1:size(edges,1)) + size(edges,1)*(edges(:,2)'-1)) = -1;
    
    % Return boundary edges.
    extEdges = zeros(size(bndry1,1),1);
    for e = 1:length(extEdges)
        extEdges(e) = find( (d0(:,bndry1(e,1)) ~= 0) & (d0(:,bndry1(e,2)) ~= 0) );
    end
    bndry1 = extEdges;
    intEdges = 1:size(d0,1);
    intEdges(extEdges) = [];
    
    ind = 2*ones(size(d0,1),1);
    ind(bndry1) = 1;
    ind = cumsum([1;ind]);
    ind = ind(1:end-1);
    
    d1 = zeros(size(Mesh.ConnectivityList,1),size(edges,1));
    triEdge = Mesh.edgeAttachments(edges);
    flatTri = [triEdge{:}];
    f1 = flatTri(ind);
    f2 = flatTri(ind(intEdges) + 1);
    pos = prod(Mesh.ConnectivityList(f1,[1,2]) == edges,2) | ...
          prod(Mesh.ConnectivityList(f1,[2,3]) == edges,2) | ...
          prod(Mesh.ConnectivityList(f1,[3,1]) == edges,2) ;
    pos = 2*pos - 1;

    d1(f1 + size(d1,1)*(0:(length(edges)-1))) = pos;
    d1(f2 + size(d1,1)*(intEdges-1)) = -pos(intEdges);

    bndry2 = find(sum(abs(d1(:,bndry1)),2) > 0);
    
end

