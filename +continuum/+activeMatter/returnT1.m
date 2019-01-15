function [ Tri, t1Edges ] = returnT1( Mesh, d1 )
    %RETURN T1. Checks the triangulation for positivity condition for all
    %internal edges. Any violators will be `flipped' provided the
    %alternative fixes the problem.
    
    int1 = 1:size(d1,2);    
    intAngles = Mesh.Angles;
    
    edge = Mesh.edges;
    triEdge = Mesh.edgeAttachments(edge);
    
    edge = edge(int1,:);
    triEdge = vertcat(triEdge{int1});
    
    faceVIndx = zeros(length(int1),2);
    [~,faceVIndx(:,1)] = max(bsxfun(@ne,Mesh.Triangulation(triEdge(:,1),:),edge(:,1)) & ...
                     bsxfun(@ne,Mesh.Triangulation(triEdge(:,1),:),edge(:,2)),[],2);
    [~,faceVIndx(:,2)] = max(bsxfun(@ne,Mesh.Triangulation(triEdge(:,2),:),edge(:,1)) & ...
                     bsxfun(@ne,Mesh.Triangulation(triEdge(:,2),:),edge(:,2)),[],2);
    
    t1Criteria = .5*( cot( intAngles(triEdge(:,1) + size(intAngles,1)*(faceVIndx(:,1)-1)) ) + ...
                      cot( intAngles(triEdge(:,2) + size(intAngles,1)*(faceVIndx(:,2)-1)) ) );
                  
    t1Edges = find(t1Criteria <= 0); 
    Tri = Mesh.Triangulation;
    
    if (~isempty(t1Edges))
        % Look for edges in the list that border the same triangle.
        oldVerts = zeros(length(t1Edges),2);
        n = 1;
        triBonds = [];
        for t1E = t1Edges'
            face1 = triEdge(t1E,1);
            face2 = triEdge(t1E,2);
            vertex1 = edge(t1E,1);
            vertex2 = edge(t1E,2);
            valency1 = sum(sum(Tri==vertex1));
            valency2 = sum(sum(Tri==vertex2));
            if (~any(ismember([vertex1,vertex2],oldVerts(:))) && valency1 > 3 && valency2 > 3)
                vertexA = Tri(face1, (Mesh.Triangulation(face1,:) ~= vertex1) & (Mesh.Triangulation(face1,:) ~= vertex2)  );
                vertexB = Tri(face2, (Mesh.Triangulation(face2,:) ~= vertex1) & (Mesh.Triangulation(face2,:) ~= vertex2)  );
                Tri(face1,Tri(face1,:) == vertex2) = vertexB;
                Tri(face2,Tri(face2,:) == vertex1) = vertexA;
                oldVerts(n,:) = [vertex1,vertex2];
            elseif (valency1 == 3 || valency2 == 3)
                triBonds = [triBonds,n];
            end
            n = n + 1;
        end
        t1Edges(triBonds) = [];
        t1Edges = int1(t1Edges); % Reindex to agree with original list.
    end
    
end

