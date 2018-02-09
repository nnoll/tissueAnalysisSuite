classdef surface < handle
    %   Generates and stores a time-series of meshes created by user 
    %   supplied dynamical equations.
    
    %---------------------------------------------------------------------
    % properties
    %---------------------------------------------------------------------

    properties (SetAccess = protected)

        a_n %Power series coefficients for a(z)
        b_n %Power series coefficients for b(z)
        x_n %Power series coefficients for x(z)
        y_n %Power series coefficients for y(z)
        Da  %Power series coefficients for the derivative of a(z)
        Db  %Power series coefficients for the derivative of b(z)
        Dx  %Power series coefficients for the derivative of x(z)
        Dy  %Power series coefficients for the derivative of y(z)
        z_B %Stores [zmin,zmax] so we know where caps are.
        
        aS
        aMu
        bS
        bMu
        xS
        xMu
        yS
        yMu
        
    end
    

    %---------------------------------------------------------------------
    % public methods
    %---------------------------------------------------------------------
    
    methods
        
        %------------------------------------------------------
        % constructor
        %------------------------------------------------------
        
        function this = surface( r, n, vargin )
            % SURFACE Creates a new reference surface.
            %
            % surface(r, n): Fits and stores a reference surface from the
            % vertex positions r. Stores an nth order polynomial.
            
            if (nargin == 2)
                % Initialize the grid in Z.
                numBins = round(size(r,1)/20);
                zG = linspace(min(r(:,3)),max(r(:,3)),numBins+1);

                a = zeros(numBins,1);
                b = zeros(numBins,1);
                x0 = zeros(numBins,1);
                y0 = zeros(numBins,1);
                
                z0 = (zG(1:(end-1)) + zG(2:end))'/2;
                badBins = [];

                for ii = 1:numBins
                    verts = (r(:,3) <= zG(ii+1)) .* (r(:,3) >= zG(ii)) == 1;
                    if (sum(verts) >= 5)
                        [ ellipse ] = continuum.surface.fitEllipse(r(verts,1),r(verts,2));
                        a(ii) = ellipse.a;
                        b(ii) = ellipse.b;
                        x0(ii) = ellipse.X0;
                        y0(ii) = ellipse.Y0;
                    else
                        badBins = horzcat(badBins,ii);
                    end
                end

                z0(badBins) = [];
                a(badBins) = [];
                b(badBins) = [];
                x0(badBins) = [];
                y0(badBins) = [];

                % Make max and min the poles.
                z0E = vertcat(zG(1),z0,zG(end));
                a = vertcat(0,a,0);
                b = vertcat(0,b,0);
%                 x0 = vertcat(0,x0,0);
%                 y0 = vertcat(0,y0,0);

                %Intialize the structure's properties.
                [this.a_n,this.aS,this.aMu] = polyfit(z0E,a,n);
                
                [this.b_n,this.bS,this.bMu] = polyfit(z0E,b,n);
                
                [this.x_n,this.xS,this.xMu] = polyfit(z0,x0,n-2);
                
                [this.y_n,this.yS,this.yMu] = polyfit(z0,y0,n-2);
                
                this.Da = polyder(this.a_n);
                this.Db = polyder(this.b_n);
                this.Dx = polyder(this.x_n);
                this.Dy = polyder(this.y_n);
                this.z_B = [zG(1),zG(end)];
                
            elseif (nargin == 3)
                this.a_n = r;
                this.b_n = n;
                this.Da = vargin{1};
                this.Db = vargin{2};
                this.z_B = vargin{3};
            end
               
        end
        
        %------------------------------------------------------
        % Evaluates the X and Y positions given Z and PHI surface coords.
        %------------------------------------------------------
        
        [X,Y] = computeCartesian( this, Z, P );  
        
        %------------------------------------------------------
        % Computes normal vector given either Z and PHI surface coords.
        % OR vertex positions (assuming they are on the surface).
        %------------------------------------------------------
        
        [ N ] = computeNormal( this, Z, P );
        
        %------------------------------------------------------
        % Computes the normalized surface basis vectors.
        %------------------------------------------------------
        
        [ e_phi, e_z ] = computeSurfaceBasis( this, Z, P );
        
        %------------------------------------------------------
        % Computes the surface coordinates.
        %------------------------------------------------------
        
        [ z, phi ] = computeSurfaceCoords( this, v );
        
        %------------------------------------------------------
        % Computes the distance to the surface for a given set of points.
        %------------------------------------------------------
        
        [ d ] = computeDistance( this, r );
        
        %------------------------------------------------------
        % Plots the fitted surface as a mesh object.
        %------------------------------------------------------
        
        plotSurface( this );
        
        %------------------------------------------------------
        % Save object.
        %------------------------------------------------------
        
%         function [saveObj] = saveobj(obj)
%             saveObj = struct('a_n',obj.a_n,'b_n',obj.b_n,'Da',obj.Da,'Db',obj.Db,'z_B',obj.z_B);
%         end
        
    end
    
%     methods (Static)
%         
%         %------------------------------------------------------
%         % Load object.
%         %------------------------------------------------------
%         
%         function [loadTopo] = loadobj(obj)
%             vargin{1} = obj.Da;
%             vargin{2} = obj.Db;
%             vargin{3} = obj.z_B;
%             loadTopo = surface.surface(obj.a_n,obj.b_n,vargin);
%         end
%         
%     end
    
end

