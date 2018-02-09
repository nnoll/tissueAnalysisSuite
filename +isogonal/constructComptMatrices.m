function [ compt, B1, B2 ] = constructComptMatrices( tri, bulkCells, extCells, bulkVerts, extVerts, Struct )
    % CONSTRUCT COMPT MATRICES 

    iCells = [bulkCells,extCells];
    iVerts = [bulkVerts,extVerts];
    
    compt = zeros(length(bulkCells),3*length(bulkVerts));
    nV = size(tri,1);
    
    for ii = 1:length(bulkCells)
        
        % Find all verts connected to cell c 
        verts1 = find(tri(:,1) == ii);
        verts2 = find(tri(:,2) == ii);
        verts3 = find(tri(:,3) == ii);
        
        compt(ii,verts1 + nV) = 1;
        compt(ii,verts1 + 2*nV) = -1;
        
        compt(ii,verts2 + 2*nV) = 1;
        compt(ii,verts2) = -1;
        
        compt(ii,verts3) = 1;
        compt(ii,verts3 + nV) = -1;
        
    end
    
    B1 = zeros(3*length(bulkVerts),length(iVerts));
    B2 = zeros(3*length(bulkVerts),length(iVerts));
    
    for t = 1:size(tri,1)
        
        v = bulkVerts(t);
        % Vertex Neighbor between 1/2
        nv12 = Struct.Cdat(iCells(tri(t,1))).nverts(ismember(...
              Struct.Cdat(iCells(tri(t,1))).nverts, ...
              Struct.Cdat(iCells(tri(t,2))).nverts));
        nv12 = nv12(nv12~=v);
        
        % Vertex Neighbor between 2/3
        nv23 = Struct.Cdat(iCells(tri(t,2))).nverts(ismember(...
              Struct.Cdat(iCells(tri(t,2))).nverts, ...
              Struct.Cdat(iCells(tri(t,3))).nverts));
        nv23 = nv23(nv23~=v);

        % Vertex Neighbor between 3/1
        nv31 = Struct.Cdat(iCells(tri(t,3))).nverts(ismember(...
              Struct.Cdat(iCells(tri(t,3))).nverts, ...
              Struct.Cdat(iCells(tri(t,1))).nverts));
        nv31 = nv31(nv31~=v);
        
        % Find indices in iVerts vector.
%         v
%         iCells(tri(t,:))
%         nv12
%         nv23
%         nv31
        v = iVerts == v;
        nv12 = iVerts == nv12(1);
        nv23 = iVerts == nv23(1);
        nv31 = iVerts == nv31(1);
        
        % Positive Bonds
        B1(t,nv31) = 1;
        B1(t,v) = -1;
        
        B1(t+nV,nv12) = 1;
        B1(t+nV,v) = -1;
        
        B1(t+2*nV,nv23) = 1;
        B1(t+2*nV,v) = -1;
        
        % Negative Bonds
        B2(t,nv12) = 1;
        B2(t,v) = -1;
        
        B2(t+nV,nv23) = 1;
        B2(t+nV,v) = -1;
        
        B2(t+2*nV,nv31) = 1;
        B2(t+2*nV,v) = -1;
    end
    
    compt = sparse(compt);
    B1 = sparse(B1);
    B2 = sparse(B2);
    
end

