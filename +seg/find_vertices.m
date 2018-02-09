function [ Vdat ] = find_vertices( L, threeFold )

    %Find position of all vertices
    if (nargin == 1 || threeFold == 0)
        CC = bwconncomp(seg.findBranchPoints( L==0 ), 8);

        S = regionprops(CC,'Centroid');
        Vertices = cell(length(S),1);
        for ii = 1:length(S)
            Vertices{ii} = round(S(ii).Centroid);
        end
    else
%         rgb(:,:,1) = (L==0);
%         rgb(:,:,2) = seg.findBranchPoints(L==0);
%         rgb(:,:,3) = 0;
%         imshow(uint8(255*rgb))
%         pause
        [y,x] = find(seg.findBranchPoints(L==0));
        Vertices = cell(length(x),1);
        for ii = 1:length(y)
            Vertices{ii} = [x(ii),y(ii)];
            CC.PixelIdxList{ii} = y(ii) + size(L,1)*(x(ii)-1);
        end
        
    end
        
    %Find neighboring cells of vertices.
    for v = 1:length(Vertices)
        nCells = 9*length(CC.PixelIdxList{v});
        for ii = 1:length(CC.PixelIdxList{v})
            nCells((9*(ii-1)+1):(9*ii)) = CC.PixelIdxList{v}(ii) + ...
            [0,-1,1,size(L,1),size(L,1)-1,size(L,1)+1,-size(L,1),-size(L,1)-1,-size(L,1)+1];
        end
        nCells = L(nCells);
        nCells = unique(nCells(nCells~=0));
        
%         if (length(nCells) == 2)
%             imshow(imdilate(L==0,strel('disk',2)))
%             hold all
%             scatter(Vertices{v}(1),Vertices{v}(2))
%             pause
%         end
            
        Vertices{v} = horzcat(uint16(Vertices{v}),uint16(nCells));
    end
    
    N = length(Vertices);
    adj = zeros(N,N);
    Vdat(N) = struct('vertxcoord',[],'vertycoord',[],'ncells',[],'nverts',[]);

    %Extract out neighboring cells.
    R = zeros(2,N);
    for ii=1:N
        Vdat(ii).vertxcoord = Vertices{ii}(1);
        Vdat(ii).vertycoord = Vertices{ii}(2);
        R(1,ii) = Vertices{ii}(1);
        R(2,ii) = Vertices{ii}(2);
        Vdat(ii).ncells = Vertices{ii}(3:end);
    end

    D = bsxfun(@plus,dot(R,R,1)',dot(R,R,1)) - 2*(R'*R);
    %Compute vertex connections via overlap of neighboring cells.
    for ii = 1:N
        for jj = ii+1:N
           if (D(ii,jj) <= 200000)
               common_mems = ismembc(Vdat(ii).ncells,Vdat(jj).ncells);
               ncommon = sum(common_mems);
               if (ncommon >= 2)
                   adj(ii,jj) = 1;
                   adj(jj,ii) = 1;
               end
           end
        end
        Vdat(ii).nverts = find(adj(ii,:)==1);
    end

    for ii = 1:length(Vdat)
        if (length(Vdat(ii).nverts) == 4)
            Vdat(ii).fourfold = 1;
        end
    end
  
end



