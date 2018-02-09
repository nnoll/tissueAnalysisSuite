function [ d0 ] = primalBonds( Tri )
    %GENERATE TRIANGULATION 
    
    % Define exterior derivatives on vertices and edges.
    edge = edges(Tri);
    d0 = zeros(size(edge,1),size(Tri.X,1));
    d0((1:size(edge,1)) + size(edge,1)*(edge(:,1)'-1)) = 1;
    d0((1:size(edge,1)) + size(edge,1)*(edge(:,2)'-1)) = -1;

end

