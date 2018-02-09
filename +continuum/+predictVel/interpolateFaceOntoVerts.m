function [ TV ] = interpolateFaceOntoVerts( Mesh, T )

    dim = size(T);
    TV = zeros([size(Mesh.X,1),dim(2:end)]);
    neighbors = Mesh.vertexAttachments;
    for v = 1:size(Mesh.X,1)
        TV(v,:) = mean(T(neighbors{v},:),1);
    end
    
end

