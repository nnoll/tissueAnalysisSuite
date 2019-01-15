function [ skeleton ] = imgSkel( Struct, X, Y, mode )
%SKELETONIZE_VERTS Takes the network topology and triangulated position of
%vertices to produce a skeletonized image.

skeleton = zeros(Y, X, length(Struct));
sed = strel('disk',1);
g = fspecial('gaussian', 4);

if (nargin == 3 || mode == 0)
    for t=1:length(Struct)
        Vlist = Struct(t).Vdat;

        adj_matrx = zeros(length(Vlist));
        for v = 1:length(Vlist)
            connverts = Vlist(v).nverts;
            kk = 1;
            for nv = connverts
                if (~adj_matrx(v,nv) && ~(ismember(1,Vlist(v).ncells) && ismember(1,Vlist(nv).ncells)) )
                    %Draw line between verts
                    [x, y] = plot.bresenham(Vlist(v).vertxcoord,Vlist(v).vertycoord,Vlist(nv).vertxcoord,Vlist(nv).vertycoord);
                    inds = y + Y*(x-1) + Y*X*(t-1);

                    skeleton(inds) = 1;
                    %Update adjacency matrix
                    adj_matrx(v,nv) = 1; adj_matrx(nv,v) = 1;
                end
                kk = kk + 1;
            end
        end
        skeleton(:,:,t) = imdilate(skeleton(:,:,t),sed);
        skeleton(:,:,t) = imfilter(skeleton(:,:,t),g);
    end
else
    for t=1:length(Struct)
        Vlist = Struct(t).Vdat;

        adj_matrx = zeros(length(Vlist));
        for v = 1:length(Vlist)
            connverts = Vlist(v).nverts;
            kk = 1;
            for nv = connverts
                if (~adj_matrx(v,nv) && (~(ismember(1,Vlist(v).ncells) || ismember(1,Vlist(nv).ncells))) )
                    %Draw line between verts
                    [x, y] = plot.bresenham(Vlist(v).vertxcoord,Vlist(v).vertycoord,Vlist(nv).vertxcoord,Vlist(nv).vertycoord);
                    inds = y + Y*(x-1) + Y*X*(t-1);

                    skeleton(inds) = 1;
                    %Update adjacency matrix
                    adj_matrx(v,nv) = 1; adj_matrx(nv,v) = 1;
                end
                kk = kk + 1;
            end
        end
        skeleton(:,:,t) = imdilate(skeleton(:,:,t),sed);
        skeleton(:,:,t) = imfilter(skeleton(:,:,t),g);
    end
    
end


end

