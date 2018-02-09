function [ n, celln ] = extract_fittedNet( T, Struct, three_verts, involved_cells )
%EXTRACT_FITTEDNET Takes in the fitted tension net and backs out the
%normalized triplet of vectors at each vertex.

%%Inputs
%1. T - Location of dual lattice vertices.
%2. Struct - Data structure.
%3. three_verts - list of vertices in the bulk.
%4. involved_cells - list of cells used for the fit.

%%Outputs
%1. n - Triplet at each vertex.
%2. celln - Ordered Triplet at each cell.

n = zeros(2,3,length(three_verts));
celln = zeros(3,length(three_verts));
R = [0, -1; 1, 0];
for ii=1:length(three_verts)
    v = three_verts(ii);
    ncells = Struct.Vdat(v).orderedncells;
    
    %Compute cell indices
    cell_ind = zeros(1,3);
    for jj=1:3
        cell_ind(jj) = find(involved_cells==ncells(jj),1);
    end
    celln(:,ii) = ncells;
     
    %Compute normalized tension edge vectors.
    n(:,1,ii) = T(:,cell_ind(2)) - T(:,cell_ind(1));
    n(:,2,ii) = T(:,cell_ind(3)) - T(:,cell_ind(2));
    n(:,3,ii) = T(:,cell_ind(1)) - T(:,cell_ind(3));
    n(:,1,ii) = n(:,1,ii)/norm(n(:,1,ii));
    n(:,2,ii) = n(:,2,ii)/norm(n(:,2,ii));
    n(:,3,ii) = n(:,3,ii)/norm(n(:,3,ii));
    
    %Rotate
    n(:,1,ii) = R*n(:,1,ii);
    n(:,2,ii) = R*n(:,2,ii);
    n(:,3,ii) = R*n(:,3,ii);
    
end

end

