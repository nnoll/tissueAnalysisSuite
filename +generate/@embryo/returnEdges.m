function [ edgeCells, edgeVerts ] = returnEdges( this, d0, d1 )
%RETURNEDGES Summary of this function goes here
%   Detailed explanation goes here

    edgeCells = zeros(size(d0,1),2);
    edgeVerts = zeros(size(d0,1),2);

    badVerts = [];
    for b = 1:size(d0,1)
        
        if (~isempty(find(d0(b,:)==1)))
            edgeCells(b,1) = find(d0(b,:)==1);
        else
            edgeCells(b,1) = 0;
        end
        
        if (~isempty(find(d0(b,:)==-1)))
            edgeCells(b,2) = find(d0(b,:)==-1);
        else
            edgeCells(b,1) =0;
        end
        
        if (sum(abs(d1(:,b))) == 2)
            edgeVerts(b,1) = find(d1(:,b)==1);
            edgeVerts(b,2) = find(d1(:,b)==-1);
        end
    end

end

