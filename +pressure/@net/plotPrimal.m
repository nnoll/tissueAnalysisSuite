function [ ] = plotPrimal( this, color, r0 )
%PLOT PRIMAL plot the segmented cells based on supplied colors
%
% r0 : probably the positions of the dual vertices or the vertices??
    
    if (nargin == 3)
        r = this.computePrimalVerts(r0);
        if (isempty(color))
            color = 'k';
        end
    else
        r = this.computePrimalVerts();
    end
    [ faces ] = this.computeFaces( r );
    
    goodFaces = sum(~isnan(faces),2) > 3;
    faces = faces(goodFaces,:);
    
    if (nargin == 1)
        patch('Vertices',r,'Faces',faces,'FaceColor','none','LineWidth',2,'EdgeColor','b')
    else
        if (ischar(color))
            patch('Vertices',r,'Faces',faces,'FaceColor','none','LineWidth',2,'EdgeColor',color)
        elseif (length(color) == 3)
            patch('Vertices',r,'Faces',faces,'FaceColor','none','LineWidth',2,'EdgeColor',color)
        else
            patch('Vertices',r,'Faces',faces,'FaceVertexCData',color(goodFaces),'LineWidth',2,'EdgeColor','b')
            shading faceted
        end
    end
    
end

