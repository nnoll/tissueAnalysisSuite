function [ L2_t, del ] = removeBadCells( Struct, L_t, mode )
%RM_BADCELLS Removes 'bubble' cells that only have two verts in the
%segmented image.

if (nargin == 2 || mode == 0)
    L2_t = L_t;
    
    for t = 1:length(Struct)
        L = L_t(:,:,t);
        L2 = L2_t(:,:,t);
        bad_cells = [];
        for c=1:length(Struct(t).Cdat)
            if (length(Struct(t).Cdat(c).nverts) == 2)
                bad_cells = horzcat(bad_cells,c);
            end
        end
        
        usedCells = [];
        for bc = bad_cells
            %Identify the cell with the closest centroid. Will absorb bubble
            %into said cell.

            ncells = Struct(t).Cdat(bc).ncells;
            d1 = sum((Struct(t).Cdat(ncells(1)).centroid.coord - Struct(t).Cdat(bc).centroid.coord ).^2);
            d2 = sum((Struct(t).Cdat(ncells(2)).centroid.coord - Struct(t).Cdat(bc).centroid.coord ).^2);
            if (d1 < d2)
                parent_cell = ncells(1);
            else
                parent_cell = ncells(2);
            end
            L2(L==bc) = parent_cell;
            tmp = L_t(:,:,t) == bc;
            tmp(L_t(:,:,t) == parent_cell) = 1;
            tmp2 = bwmorph(tmp,'bridge');
            mem = (tmp2 - tmp);

            usedCells = [bc,ncells];
            %Remove vertices from mem
            mem(Struct(t).Vdat(Struct(t).Cdat(bc).nverts(1)).vertycoord,Struct(t).Vdat(Struct(t).Cdat(bc).nverts(1)).vertxcoord) = 0;
            mem(Struct(t).Vdat(Struct(t).Cdat(bc).nverts(2)).vertycoord,Struct(t).Vdat(Struct(t).Cdat(bc).nverts(2)).vertxcoord) = 0;
            L2(mem==1) = parent_cell;

        end
        usedCells = unique(usedCells(usedCells~=1));
        %Relabel cells to account for deletions.
        nn = 0;
        for bc = bad_cells
            tmp = L2 > (bc-nn);
            L2(tmp) = L2(tmp) - 1;
            nn = nn + 1;
        end
        
        Lbw = L2 == 0;
        Lbw2 = bwmorph(bwmorph(Lbw,'diag'),'fill');
%         Lbw2(imdilate(L==1,strel('disk',3))) = 0;
        newVerts = find((1-Lbw).*bwmorph(Lbw2,'skel'));
%         newVerts = bwmorph(bwmorph(Lbw,'skel'),'endPoints');
%         L(newVerts)
%         usedCells
        newVerts(~ismember(L(newVerts),usedCells)) = [];
%         Lbw2 = zeros(size(Lbw2));
%         Lbw2(newVerts ) = 1;
%         newVerts = Lbw2;
%         rgb(:,:,1) = Lbw;
%         rgb(:,:,2) = newVerts;
%         rgb(:,:,3) = imdilate(newVerts,strel('disk',10));
%         imshow(uint8(255*rgb))
%         pause
        if (~isempty(bad_cells))
            L2(newVerts) = 0;
        end
        
        L2_t(:,:,t) = L2;
        
    end
else
    L2 = L_t;
    del = 0;
    bad_cells = [];
    
    border = zeros(size(L_t));
    border(1,:) = 1;
    border(size(border,1),:) = 1;
    border(:,1) = 1;
    border(:,size(border,2)) = 1;
    
    for c = 1:length(Struct.Cdat)
        if (length(Struct.Cdat(c).nverts) == 2)
            % Make sure its not a border cell
            if ( sum(sum( (L_t == c) .* border)) == 0 )
                bad_cells = horzcat(bad_cells,c);
            end
        end
    end
    
    if (~isempty(bad_cells))
        del = 1;
    end
    
    for bc = bad_cells
        %Identify the cell with the closest centroid. Will absorb bubble
        %into said cell.

        ncells = Struct.Cdat(bc).ncells;
        d1 = sum((Struct.Cdat(ncells(1)).centroid.coord - Struct.Cdat(bc).centroid.coord ).^2);
        d2 = sum((Struct.Cdat(ncells(2)).centroid.coord - Struct.Cdat(bc).centroid.coord ).^2);
        
        if (d1 < d2)
            parent_cell = ncells(1);
        else
            parent_cell = ncells(2);
        end
        
        L2(L_t==bc) = parent_cell;
        tmp = L_t == bc;
        tmp(L_t == parent_cell) = 1;
        tmp2 = bwmorph(tmp,'bridge');
        mem = (tmp2 - tmp);

        %Remove vertices from mem
        mem(Struct.Vdat(Struct.Cdat(bc).nverts(1)).vertycoord,Struct.Vdat(Struct.Cdat(bc).nverts(1)).vertxcoord) = 0;
        mem(Struct.Vdat(Struct.Cdat(bc).nverts(2)).vertycoord,Struct.Vdat(Struct.Cdat(bc).nverts(2)).vertxcoord) = 0;
        L2(mem==1) = parent_cell;

    end

    %Relabel cells to account for deletions.
    nn = 0;
    for bc = bad_cells
        tmp = L2 > (bc-nn);
        L2(tmp) = L2(tmp) - 1;
        nn = nn + 1;
    end

    L2_t = L2;
    
end

end

