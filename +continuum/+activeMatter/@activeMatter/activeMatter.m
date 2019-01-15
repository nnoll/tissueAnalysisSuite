classdef activeMatter < handle
    %   Generate a primal structure to be used for discrete simulation of 
    %   continuum ATN model.
    %
    
    %---------------------------------------------------------------------
    % Properties
    %---------------------------------------------------------------------

    properties (SetAccess = protected)

        Mesh    % FELICITY Mesh structure.
        surf    % Object containing all relevant information about the surface of constraint.
        
        d0      % Exterior derivative on 0-forms.
        d1      % Exterior derivative on 1-forms.
        Lp      % Length of edges in the primal triangulation.
        Ld      % Length of edges in the dual CD.
        Ap      % Area of triangles in the primal triangulation.
        Ad      % Area of cells in the dual CD.
        
        e1      % First basis vector for each vertex.
        e2      % Second basis vector for each vertex.
        K       % Curvature vector.  
        kappa   % Gauss Curvature.
                
    end


    %---------------------------------------------------------------------
    % Public methods
    %---------------------------------------------------------------------
    
    methods
        
        %------------------------------------------------------
        % Constructor
        %------------------------------------------------------
        
        function this = activeMatter( arg1, arg2, arg3 )
            % activeMatter Creates a new primal 
            %
            % activeMatter(). Initializes default properties - hexagonal lattice 
            % in the plane.
            % activeMatter(N). N controls the number of points per linear row in
            % the plane.
            
            if (nargin == 0)
                [v,Tri] = read_ply('./embryo.ply');
%                 v = v + .0075*randn(size(v));
%                 [az,el] = cart2sph(v(:,1),v(:,2),v(:,3));
%                 el = pi/2 - el;
%                 v = [sin(el).*cos(az),sin(el).*sin(az),cos(el)];
%                 [Tri] = convhull(v(:,1),v(:,2),v(:,3));
%                 v(:,3) = -v(:,3);
                this.surf = continuum.surface.surface(v,8);
            elseif (nargin == 2)
                v = arg1;
                Tri = arg2;
                this.surf = surface.surface(v,8);
                [ z, phi ] = this.surf.computeSurfaceCoords( v );
                [x,y] = this.surf.computeCartesian( z, phi );
                v = [x,y,v(:,3)];
            else                
                % Need to project v onto the surface.
                v = arg1;
                [ z, phi ] = arg3.computeSurfaceCoords( v );
                [x,y] = arg3.computeCartesian( z, phi );
                v = [x,y,v(:,3)];
                
                Tri = arg2;
                this.surf = arg3;
            end
            
            this.Mesh = MeshTriangle(Tri,v,'Primal');
            
            [ this.d0, this.d1 ] = continuum.activeMatter.generatePrimal( this.Mesh );
            [ this.Lp, this.Ld, this.Ap, this.Ad, this.K, this.kappa ] = continuum.activeMatter.computeMetricQuantities( this.Mesh, sparse(this.d0), sparse(this.d1) );
%             [az,el] = cart2sph(this.Mesh.X(:,1),this.Mesh.X(:,2),this.Mesh.X(:,3));
%             el = pi/2-el;

%             this.e2 = [-sin(az),cos(az),zeros(size(az))];
%             this.e1 = [cos(el).*cos(az),cos(el).*sin(az),-sin(el)];
            [ Z, Phi ] = this.surf.computeSurfaceCoords( v );
            [ this.e2, this.e1 ] = this.surf.computeSurfaceBasis( Z, Phi );
            
            % Orthonormalize basis vectors
            this.e1 = this.e1 - bsxfun(@times,dot(this.e1,this.e2,2),this.e2);
            this.e1 = bsxfun(@rdivide,this.e1,sqrt(sum(this.e1.^2,2)));

            % Remove all components of e1 and e2 that lie along the
            % curvature direction.
            n = bsxfun(@rdivide,this.K,sqrt(sum(this.K.^2,2)));
            this.e1 = this.e1 - bsxfun(@times,dot(this.e1,n,2),n);
            this.e2 = this.e2 - bsxfun(@times,dot(this.e2,n,2),n);
  
            this.e1 = bsxfun(@rdivide,this.e1,sqrt(sum(this.e1.^2,2)));
            this.e2 = this.e2 - bsxfun(@times,dot(this.e2,this.e1,2),this.e1);
            this.e2 = bsxfun(@rdivide,this.e2,sqrt(sum(this.e2.^2,2)));

        end
        
    end
               
end

