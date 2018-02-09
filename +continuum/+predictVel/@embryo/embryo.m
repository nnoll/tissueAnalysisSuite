classdef embryo < handle
    %   Generate a primal structure to be used for discrete simulation of 
    %   continuum ATN model.
    %
    
    %---------------------------------------------------------------------
    % Properties
    %---------------------------------------------------------------------

    properties (SetAccess = protected)

        emb % Static mesh.
        Sm  % Smoothing kernel.
        
        Div % Divergence of Tensor

        v % Observed velocity field.
                
    end


    %---------------------------------------------------------------------
    % Public methods
    %---------------------------------------------------------------------
    
    methods
        
        %------------------------------------------------------
        % Constructor
        %------------------------------------------------------
        
        function this = embryo( )
            % embryo Creates a new embryo to work on. 

            this.emb = continuum.activeMatter.activeMatter( );
            this.emb = this.emb.smoothMesh(10);
            
            this.Sm = -.25*this.emb.d0'*diag(this.emb.Ld./this.emb.Lp)*this.emb.d0;
            this.Sm = expm(this.Sm)^5;
            
        end
        
    end
               
end


