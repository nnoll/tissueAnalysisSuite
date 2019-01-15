function [ ] = skel( Struct, color, mode )
%SKELETONIZE_VERTS Takes the network topology and triangulated position of
%vertices to produce a skeletonized image.

    if (nargin == 1)
        color = 'r';
    end
    
    Vlist = Struct.Vdat;

    adj_matrx = zeros(length(Vlist));
    hold all
    for v = 1:length(Vlist)
        
        connverts = Vlist(v).nverts;
        kk = 1;
        if (mode == 1)
            for nv = connverts
                if (~adj_matrx(v,nv) && ~(ismember(1,Vlist(v).ncells) && ismember(1,Vlist(nv).ncells)) )
                    %Draw line between verts
                    plot([Vlist(v).vertxcoord,Vlist(nv).vertxcoord],[Vlist(v).vertycoord,Vlist(nv).vertycoord],'Color',color,'LineWidth',2);

                    %Update adjacency matrix
                    adj_matrx(v,nv) = 1; adj_matrx(nv,v) = 1;
                end
                kk = kk + 1;
            end
        else
            for nv = connverts
                if ( ~adj_matrx(v,nv) )
                    %Draw line between verts
                    plot([double(Vlist(v).vertxcoord),double(Vlist(nv).vertxcoord)],...
                         [double(Vlist(v).vertycoord),double(Vlist(nv).vertycoord)],'Color',color,'LineWidth',2);

                    %Update adjacency matrix
                    adj_matrx(v,nv) = 1; adj_matrx(nv,v) = 1;
                end
                kk = kk + 1;
            end
        end
    end
    
end

