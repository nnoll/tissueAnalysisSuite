function [ ] = plotDual( this, mode )
    %PLOTDUAL 

    Xd = this.Mesh.circumcenters;
    
    if (nargin == 1 || strcmp(mode,'2D'))
        [ z, ph ] = this.surf.computeSurfaceCoords( Xd );
        z = (z - min(z))/(max(z) - min(z));
        ph(ph<0) = - pi - ph(ph<0);
        ph(ph>0) = pi - ph(ph>0);
        
        for e = 1:size(this.d1,2)
            dualVerts = (this.d1(:,e) ~= 0);
            if (abs(diff(ph(dualVerts))) < pi )
                line(z(dualVerts),ph(dualVerts),'Color',[0,0,1],'LineWidth',1);
            end
        end
        axis([0,1,-pi,pi])
    else
        for e = 1:size(this.d1,2)
            dualVerts = (this.d1(:,e) ~= 0);
            line(Xd(dualVerts,1),Xd(dualVerts,2),Xd(dualVerts,3),'Color',[0,0,1],'LineWidth',1);
        end
    end
    
end

