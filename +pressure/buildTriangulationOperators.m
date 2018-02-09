function [ d0, d1, bndry0, bndry1, bndry2 ] = buildTriangulationOperators( q, Tri )
    % BUILD TRIANGULATION OPERATORS
    
    Mesh = triangulation(Tri,q);
    [ d0, d1, bndry0, bndry1, bndry2 ] = generate.primal( Mesh );
    
end

