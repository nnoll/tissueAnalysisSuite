classdef net < handle
    %   Partitions bonds into unique identities to easily access their
    %   properties such as stored chemistry and lengths, etc.
    
    %---------------------------------------------------------------------
    % Properties
    %---------------------------------------------------------------------

    properties (SetAccess = protected)

        q     % 2D Positions of Triangulation
        theta % Z Positions of Triangulation
        p     % Pressure 
        
        tri   % Network Connectivity. 
        pTri  % Network Connectivity so that pressure weighted vertices are ordered CCW.
        d0    % Exterior derivative on verts.
        d1    % Exterior derivative on edges.
        
        cellLabels % Useful for data. Corresponds to labels in data struct.
        vertexLabels
    end


    %---------------------------------------------------------------------
    % Public methods
    %---------------------------------------------------------------------
    
    methods
        
        %------------------------------------------------------
        % Constructor
        %------------------------------------------------------
        
        function this = net( q, Theta, p, Tri, cellLabels, vertexLabels )
            
            this.q = q;
            this.theta = Theta;
            this.p = p;
            
            if (nargin >= 5)
                this.cellLabels = cellLabels;
            end
            
            if (nargin == 6)
                this.vertexLabels = vertexLabels;
            end
%             this.tri = Tri;

            pQ = bsxfun(@times,this.p,this.q);
            this.pTri = zeros(size(Tri));
            for t = 1:size(Tri,1)
                deltaR = pQ(Tri(t,:),:);
                deltaR = bsxfun(@minus,deltaR,mean(deltaR,1));
                angle = atan2(deltaR(:,2),deltaR(:,1));
                angle(angle<0) = angle(angle<0) + 2*pi;
                [~,ind] = sort(angle);
                this.pTri(t,:) = Tri(t,ind);
                % Make sure min index is first.
                [~,ind] = min(this.pTri(t,:));
                this.pTri(t,:) = circshift(this.pTri(t,:),[0,1-ind]);
            end
            
            this.tri = zeros(size(Tri));
            for t = 1:size(Tri,1)
                deltaR = this.q(Tri(t,:),:);
                deltaR = bsxfun(@minus,deltaR,mean(deltaR,1));
                angle = atan2(deltaR(:,2),deltaR(:,1));
                angle = mod(angle,2*pi);
                [~,ind] = sort(angle);
                this.tri(t,:) = Tri(t,ind);
                % Make sure min index is first.
                [~,ind] = min(this.tri(t,:));
                this.tri(t,:) = circshift(this.tri(t,:),[0,1-ind]);
            end
            
            [ this.d0, this.d1 ] = pressure.buildTriangulationOperators( this.q, this.tri );
                        
        end
        
    end
               
end

