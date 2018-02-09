classdef embryo < handle
    %   Generate a primal structure to be used for discrete simulation of 
    %   continuum ATN model.
    %
    
    %---------------------------------------------------------------------
    % Properties
    %---------------------------------------------------------------------

    properties (SetAccess = protected)

        Mesh    % FELICITY Mesh structure.
        d0      % Exterior derivative on 0-forms.
        d1      % Exterior derivative on 1-forms.
        Lp      % Length of edges in the primal triangulation.
        Ld      % Length of edges in the dual CD.
        Ap      % Area of triangles in the primal triangulation.
        Ad      % Area of cells in the dual CD.
        
        e1      % First basis vector for each vertex.
        e2      % Second basis vector for each vertex.
        K       % Curvature vector.s  
  
    end


    %---------------------------------------------------------------------
    % Public methods
    %---------------------------------------------------------------------
    
    methods
        
        %------------------------------------------------------
        % Constructor
        %------------------------------------------------------
        
        function this = embryo( v, Tri )
            % activeMatter Creates a new primal 
            %
            % activeMatter(). Initializes default properties - hexagonal lattice 
            % in the plane.
            % activeMatter(N). N controls the number of points per linear row in
            % the plane.
            %
            
            this.Mesh = triangulation(Tri,v);
            
            [ this.d0, this.d1 ] = generate.primal( this.Mesh );
            [ this.Lp, this.Ld, this.Ap, this.Ad, this.K ] = this.computeMetricQuantities();
            
            Phi = atan2(v(:,2),v(:,1));
            Theta = atan2(sqrt(v(:,1).^2+v(:,2).^2),v(:,3));
            
            this.e1 = [cos(Phi).*cos(Theta),sin(Phi).*cos(Theta),-sin(Theta)];
            this.e2 = [-sin(Phi),cos(Phi),zeros(size(v,1),1)];
            
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

