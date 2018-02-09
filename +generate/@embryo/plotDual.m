function [ ] = plotDual( this, mode, color )
    %PLOTDUAL 

    Xd = this.Mesh.circumcenter;

    if (nargin < 3)
        color = [0,0,1];
    end
    
    if (nargin == 1 || strcmp(mode,'2D'))
        [Xd(:,1),Xd(:,2)] = cart2sph(Xd(:,1),Xd(:,2),Xd(:,3));
        Xd(:,2) = pi/2 - Xd(:,2);
        for e = 1:size(this.d1,2)
            dualVerts = (this.d1(:,e) ~= 0);
            if (abs(diff(Xd(dualVerts,1))) < pi && abs(diff(Xd(dualVerts,2))) < pi/2)
                line(Xd(dualVerts,1),Xd(dualVerts,2),'Color',color,'LineWidth',1);
            end
        end
        axis([-pi,pi,0,pi])
    else
        for e = 1:size(this.d1,2)
            dualVerts = (this.d1(:,e) ~= 0);
            line(Xd(dualVerts,1),Xd(dualVerts,2),Xd(dualVerts,3),'Color',color,'LineWidth',1);
        end
    end
    
end

