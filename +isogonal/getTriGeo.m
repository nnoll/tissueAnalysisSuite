function [ Q, S, dC, dV, bulkVerts, bulkCells ] = getTriGeo( Struct )
%GET_STRUCTUREFACTORS Will take in a data structure of cells and verts and
%fit to a tension graph. Returns the angles and bond triplets.

    %Inputs
    %1. Struct - Cell and vertex data
    %2. L - watershedded image

    %Outputs
    %1. T - Normalized tension on every bond.
    %2. S - Normalized area dual to every bulk vertex.
    %3. dC - diff operator on cells.
    %4. dV - diff operator on verts.
    %5. bulkVerts - vertices considered in graph.
    %6. bulkCells - cells considered in graph.
    
    [ Q, ~, dC, dV, bulkVerts, bulkCells ] = isogonal.fitTensionGraph( Struct );
    
    S = zeros(length(bulkVerts),1);
    for v = 1:length(bulkVerts)
       edges = (dV(:,v) ~= 0); 
       verts = Q(sum(abs(dC(edges,:))) > 0,:);
       S(v) = polyarea(verts(:,1),verts(:,2));
    end

end

