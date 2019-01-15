function [ ] = plotTri( this, color )
    % PLOT TRI

    if (nargin == 1)
        color = 'r';
    end
    patch('Vertices',this.q,'Faces',this.tri,'FaceColor','none','EdgeColor',color,'LineWidth',2)
    
end

