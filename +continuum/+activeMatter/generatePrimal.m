function [ d0, d1 ] = generatePrimal( Mesh )
    %GENERATE TRIANGULATION 
    
    % Define exterior derivatives on vertices and edges.
    edges = Mesh.edges;
    d0 = zeros(size(edges,1),size(Mesh.X,1));
    d0((1:size(edges,1)) + size(edges,1)*(edges(:,1)'-1)) = 1;
    d0((1:size(edges,1)) + size(edges,1)*(edges(:,2)'-1)) = -1;
    
    intEdges = 1:size(d0,1);
    
    ind = 2*ones(size(d0,1),1);
    ind = cumsum([1;ind]);
    ind = ind(1:end-1);
    
    d1 = zeros(size(Mesh.Triangulation,1),size(edges,1));
    triEdge = Mesh.edgeAttachments(edges);
    flatTri = [triEdge{:}];
    f1 = flatTri(ind);
    f2 = flatTri(ind(intEdges) + 1);
    pos = prod(Mesh.Triangulation(f1,[1,2]) == edges,2) | ...
          prod(Mesh.Triangulation(f1,[2,3]) == edges,2) | ...
          prod(Mesh.Triangulation(f1,[3,1]) == edges,2) ;
    pos = 2*pos - 1;

    d1(f1 + size(d1,1)*(0:(length(edges)-1))) = pos;
    d1(f2 + size(d1,1)*(intEdges-1)) = -pos(intEdges);
    
end

