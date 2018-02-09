function [ Struct ] = findInteriorHoles( Struct, L )
    % FIND INTERIOR HOLES 

    for t = 1:length(Struct)
        if (nargin == 2)
            S = regionprops(L(:,:,t),'Area');
            A = [S.Area];
            for c = 1:length(Struct(t).Cdat)
                verts = Struct(t).Cdat(c).nverts; % Not ordered as that assumed convexity!
                if (length(verts) > 2 && c > 1)
                    % Compute all interior angles
                    r = [Struct(t).Vdat(verts).vertxcoord;Struct(t).Vdat(verts).vertycoord]';
                    K = convhull(r(:,1),r(:,2));
                    if ((length(K)-1) < length(verts) && A(c) > prctile(A,75))
%                         imshow(L(:,:,t)==0)
%                         hold on
%                         scatter(r(:,1),r(:,2),'b','filled')
%                         scatter(r(K,1),r(K,2),'ro','LineWidth',2)
%                         pause
                        Struct(t).Cdat(c).hole = 1;
                    else
                        Struct(t).Cdat(c).hole = 0;
                    end
                else
                    Struct(t).Cdat(c).hole = 1;
                end
            end
        else
            Struct(t).Cdat(1).hole = 1;
            for c = 2:length(Struct(t).Cdat)
                Struct(t).Cdat(c).hole = 1-Struct(t).Cdat(c).all_threefold;
            end
        end
    end
    
end

