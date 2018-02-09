function [ Struct ] = findInteriorHoles( Struct, L )
    % FIND INTERIOR HOLES 

    for t = 1:length(Struct)
        S = regionprops(L(:,:,t),'Area');
        A = [S.Area];
        for c = 2:length(Struct(t).Cdat)
            verts = Struct(t).Cdat(c).nverts; % Not ordered as that assumed convexity!
            if (length(verts) > 1)
                orderedVerts = zeros(size(verts));
                orderedVerts(1) = verts(1);
                for ii = 2:length(verts)
                    nverts = Struct(t).Vdat(orderedVerts(ii-1)).nverts;
                    cverts = nverts(ismember(nverts,verts));
                    cverts = cverts(~ismember(cverts,orderedVerts));
                    if (~isempty(cverts))
                        orderedVerts(ii) = cverts(1);
                    else
                        break
                    end
                end

                if (sum(orderedVerts==0) == 0)
                    % Compute all interior angles
                    r = [Struct(t).Vdat(orderedVerts).vertxcoord;Struct(t).Vdat(orderedVerts).vertycoord]';
                    rp = [Struct(t).Vdat(circshift(orderedVerts,[0,-1])).vertxcoord;Struct(t).Vdat(circshift(orderedVerts,[0,-1])).vertycoord]';
                    rm = [Struct(t).Vdat(circshift(orderedVerts,[0,1])).vertxcoord;Struct(t).Vdat(circshift(orderedVerts,[0,1])).vertycoord]';

                    rB1 = rp - r;
                    rB2 = rm - r;

%                     theta = atan2(rB1(:,1).*rB2(:,2) - rB1(:,2).*rB2(:,1), dot(rB1,rB2,2));
                    theta = mod(atan2(rB1(:,2),rB1(:,1))',2*pi);
                    theta = 2*pi - mod(pi - diff([theta(end),theta]),2*pi);
                    if (c == 915)
                        scatter(r(:,1),r(:,2),'ro')
                        hold on
                        scatter(r(1,1),r(1,2),'g','filled')
                        scatter(r(2,1),r(2,2),'c','filled')
                        plot([r(:,1),rp(:,1)],[r(:,2),rp(:,2)],'b')
                        pause
                        rB1
                        rB2
                        180*theta/pi
                        pause
                    end
                    if (any(theta > pi) && A(c) > prctile(A,75))
                        Struct(t).Cdat(c).hole = 1;
                    else
                        Struct(t).Cdat(c).hole = 0;
                    end
                else
                    Struct(t).Cdat(c).hole = 1;
                end
            else
                Struct(t).Cdat(c).hole = 1;
            end
        end
    end
    
end

