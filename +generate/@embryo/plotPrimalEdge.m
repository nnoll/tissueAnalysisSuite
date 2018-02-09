function [ ] = plotPrimalEdge( this, edgeColor, norm, mode )
    % PLOT PRIMAL EDGE

    colormap = jet(512);
    linIdx = linspace(0,1,512);
    if (nargin == 2)
        norm = 1;
    end
    if (nargin <= 3)
        mode = '3D';
    end
    
    if (norm)
        if (std(edgeColor) > 0)
            edgeColor = (edgeColor - min(edgeColor))/(max(edgeColor)-min(edgeColor));
        else
            edgeColor = zeros(size(edgeColor));
        end
    end
    
    desiredColor = zeros(length(edgeColor),3);
    for ii = 1:3
        desiredColor(:,ii) = interp1(linIdx,colormap(:,ii),edgeColor);
    end

    hold on 
    if (strcmp(mode,'3D'))
        for e = 1:size(this.d0,1)
            verts = (this.d0(e,:) ~= 0);
            line(this.Mesh.X(verts,1),this.Mesh.X(verts,2),this.Mesh.X(verts,3),'Color',desiredColor(e,:),'LineWidth',1);
        end
    else
        [ph,th] = cart2sph(this.Mesh.X(:,1),this.Mesh.X(:,2),this.Mesh.X(:,3));
        ph(ph<0) = - pi - ph(ph<0);
        ph(ph>0) = pi - ph(ph>0);
        th = pi/2-th;
        for e = 1:size(this.d0,1)
            verts = (this.d0(e,:) ~= 0);
            if (abs(diff(ph(verts))) < pi/2)
                line(ph(verts),th(verts),'Color',desiredColor(e,:),'LineWidth',1);
            end
            set(gca,'XLim',1.05*[-pi,pi],'YLim',1.05*[0,pi]);
        end
    end
    
    colorbar()

end

