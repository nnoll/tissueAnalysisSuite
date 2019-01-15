function [ ] = plotPrimal( this, color, mode )
    % PLOT PRIMAL 
    
    if (nargin == 1 || ~strcmp(mode,'2D'))
        if (nargin == 1 || isempty(color))
            patch('Faces',this.Mesh.Triangulation,'Vertices',this.Mesh.X,'FaceColor','w','EdgeColor',[0,0,0]);
        else
            patch('Faces',this.Mesh.Triangulation,'Vertices',[this.Mesh.X(:,3),this.Mesh.X(:,2),-this.Mesh.X(:,1)],'FaceVertexCData',color,'EdgeColor',[0,0,0]);
            shading faceted
        end
        axis equal
    else
        [ z, ph ] = this.surf.computeSurfaceCoords( this.Mesh.X );
        z = 2*pi*.8482*(z - min(z))/(max(z) - min(z));
        
        ph = -ph;
        ph = mod(ph+pi+2.892,2*pi) - pi;
        
        Tri_phi = ph(this.Mesh.Triangulation); 
        Tri_phi = [Tri_phi,Tri_phi(:,1)];
        Tri_phi = diff(Tri_phi,1,2);
        goodTris = (sum(abs(Tri_phi)>pi/2,2)) < 1;
        
        if (~isempty(color))
            if (length(color) == length(z))
                patch('Faces',this.Mesh.Triangulation(goodTris,:),'Vertices',[z,ph],'FaceVertexCData',color,'EdgeColor',[0,0,0]);
                shading faceted
            else
                patch('Faces',this.Mesh.Triangulation(goodTris,:),'Vertices',[z,ph],'FaceVertexCData',color(goodTris),'EdgeColor',[0,0,0]);
            end
        else
            patch('Faces',this.Mesh.Triangulation(goodTris,:),'Vertices',[z,ph],'FaceColor','w','EdgeColor',[0,0,0]);
        end
        set(gca,'XLim',[0,2*pi*.8482],'YLim',[-pi,pi]);
    end

end

