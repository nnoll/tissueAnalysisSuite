function [Ltemp,Struct] = generate_structs(Ltrack, bond, clear_border, threefold, very_far)
    % Skeletonize the label matrix to create vertices and connectivity
    % 
    % 
    % Parameters
    % ----------
    % Ltrack : int matrix (N x M x #timepoints or zslices)
    %   segmentation of the image
    % bond : int (0 or 1)
    %   if 1, then find bonds. This is slow, so consider avoiding (set
    %   bond=0).
    % clear_border : bool or int (0/1)
    %   Typical segmentation has cells running into boundary. This changes
    %   all of the boundary cells to cell 1 (the external cell).
    % threefold : bool or int (0/1)
    %   Inverse expects threefold vertices, so this moves 4-fold vertices
    %   apart.
    % very_far : int or float
    %   How far away for two vertices to definitely not be connected.
    %   Larger than the size of a cell.
    %
    % Returns
    % -------
    % Ltemp : label matrix
    %   edges will have zero pixel value. Cells will have integer label.
    %   The ordering of cells goes down then back up (and then back down) 
    %   as we move to the right.
    % Struct : data structure
    %   skeletonization of the tissue with fields
    %       vdat: 
    %           nverts : neighbor list for vertices
    %           ncells : cells that define each vertex?
    %           vertxcoord : column of data in which vertex lives
    %           vertycoord : row of data in which each vertex lives
    %       cdat : cell data    
    %           ncells : indices of neighboring cells
    %           nverts : vertices that define the cell
    %           centroid.coord(1) x position of cell centroid
    %           centroid.coord(2) y position of cell centroid
    %       bdat : 
    %           nverts
    %           ncells
    %           pix : linear indices of the pixels associated with that 
    %               bond
    %           chem : str 
    %               which color channel the segment or intensity is in?
    %
    Ltemp=Ltrack;
    
    if (nargin == 3)
        threefold = 0;
        very_far = max(size(Ltrack(:,:,1))) ;
    elseif (nargin == 4)
        very_far = max(size(Ltrack(:,:,1))) ;
    end
    
    if (clear_border)
        for ii=1:length(Ltemp(1,1,:))
            L=Ltemp(:,:,ii);

            L1=imclearborder(L);
            L2=imsubtract(L,L1)>0;
            L3=imdilate(L2,strel('disk',1));
            L4=imerode(L3,strel('disk',1));
            L1=L1+1;L1(L1==1)=0;L1(L4>0)=1;

            L1(1:end,1)=1;L1(1:end,end)=1;L1(1,1:end)=1;L1(end,1:end,1)=1;

            %Remove spurs on exterior edge
            img = L1==0; 
            img2 = bwmorph(img,'spur');
            img3 = imsubtract(img,img2);
            L1(find(img3)) = 1;
            Ltemp(:,:,ii)=L1;

            clear L1 L2 L3 L4 
        end  
    end

    for ii=1:length(Ltemp(1,1,:))
        disp('Computing Struct for slice ', num2str(ii)])
        L = Ltemp(:,:,ii);
        Struct(ii) = seg.create_Cdat_Vdat_initial(L,bond,threefold, very_far);
        for v = 1:length(Struct(ii).Vdat)
            Struct(ii).Vdat(v).vertxcoord = double(Struct(ii).Vdat(v).vertxcoord);
            Struct(ii).Vdat(v).vertycoord = double(Struct(ii).Vdat(v).vertycoord);
        end
    end 

end
    
