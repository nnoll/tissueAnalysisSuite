function [ ] = plotVectorField( this, VF, mode, color )

    if (nargin == 2 || strcmp(mode,'3D'))
%         this.plotPrimal();
%         subplot(2,1,1)
        patch('Vertices',[this.Mesh.X(:,3),this.Mesh.X(:,2),-this.Mesh.X(:,1)], ...
              'Faces', this.Mesh.Triangulation,'FaceColor','w','EdgeColor','w');
        hold on
        axis equal
        quiver3(this.Mesh.X(:,3),this.Mesh.X(:,2),-this.Mesh.X(:,1),VF(:,3),VF(:,2),-VF(:,1),5,color);
        view([0,0]);
        
%         patch('Vertices',[this.Mesh.X(:,1),this.Mesh.X(:,2),this.Mesh.X(:,3)], ...
%               'Faces', this.Mesh.Triangulation,'FaceColor','w','EdgeColor','w');
%         hold on
%         axis equal
%         quiver3(this.Mesh.X(:,1),this.Mesh.X(:,2),this.Mesh.X(:,3),VF(:,1),VF(:,2),VF(:,3),5,color);
%         view([0,0]);
        
%         subplot(2,1,2)
%         patch('Vertices',[-this.Mesh.X(:,3),this.Mesh.X(:,2),-this.Mesh.X(:,1)], ...
%               'Faces', this.Mesh.Triangulation);
%         hold on
%         quiver3(-this.Mesh.X(:,3),this.Mesh.X(:,2),-this.Mesh.X(:,1),-VF(:,3),VF(:,2),-VF(:,1),5,color);
%         view([180,0]);

    else
        [z,phi] = this.surf.computeSurfaceCoords(this.Mesh.X);
        [ephi,ez] = this.surf.computeSurfaceBasis(z,phi); 
        z = 2*pi*.8482*(z-min(z))/(max(z)-min(z));
        
        phi = -phi;
        phi = mod(phi+pi+2.892,2*pi) - pi;
        
        phi = 2*phi;
        
        % Build metric
        g = zeros(2,2,size(ephi,1));
        g(1,1,:) = dot(ephi,ephi,2);
        g(1,2,:) = dot(ephi,ez,2);
        g(2,1,:) = dot(ez,ephi,2);
        g(2,2,:) = dot(ez,ez,2);
        
        v = [dot(VF,ephi,2)';dot(VF,ez,2)'];
        
        VF = arrayfun(@(i) g(:,:,i) \ v(:,i) ,1:size(g,3), 'UniformOutput', false);
        VF = cat(2,VF{:})';
        
%         phi(phi<0) = - pi - phi(phi<0);
%         phi(phi>0) = pi - phi(phi>0);
        
%         Sm = sparse(-.25*this.d0'*diag(this.Ld./this.Lp)*this.d0);
%         Sm = expm(Sm)^10;

%         this.plotPrimal([],'2D')
%         hold on
        quiver(2*z,phi,.8482*VF(:,2),-2*VF(:,1),1.2)
        set(gca,'XLim',2*.8482*[0,2*pi],'YLim',2*[-pi,pi])
        axis equal
    end
    
end

