function [ Tri ] = computeTriangulation( Adj )
    %COMPUTETRIANGULATION Summary of this function goes here
    %   Detailed explanation goes here

    adjSq = Adj*Adj;
    triangles = Adj .* adjSq;
    [v1,v2] = ind2sub(size(triangles),find(triangles));
    Tri = zeros(sum(triangles(:))/6,3);
    f = 1;
    for ii = 1:length(v1)
        neigh = find(Adj(v1(ii),:).*Adj(v2(ii),:));
        for n = neigh
            verts = sort([v1(ii),v2(ii),n]);
            triExist = sum( (Tri(:,1) == verts(1)) .* (Tri(:,2) == verts(2)) .* (Tri(:,3) == verts(3)) );
            if (triExist == 0)
                Tri(f,:) = verts;
                f = f + 1;
            end
        end
    end
    if (f < sum(triangles(:))/6)
        Tri(f:size(Tri,1),:) = [];
    end

end

