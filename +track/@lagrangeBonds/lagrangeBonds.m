classdef lagrangeBonds < handle
    %   Partitions bonds into unique identities to easily access their
    %   properties such as stored chemistry and lengths, etc.
    
    %---------------------------------------------------------------------
    % Properties
    %---------------------------------------------------------------------

    properties (SetAccess = protected)

        times % Times bond existed in.
        label % Label in the data structure at the given time.
        delta % Bond must be tracked longer than delta to count.
        cells % Cell IDs that bond are defined relative to
        
        r     % Stored bond length
        myo   % Stored myo concentration.
        theta % Stored bond angle
        
    end


    %---------------------------------------------------------------------
    % Public methods
    %---------------------------------------------------------------------
    
    methods
        
        %------------------------------------------------------
        % Constructor
        %------------------------------------------------------
        
        function this = lagrangeBonds( bonds, delta )

            this.delta = delta;
            badTracks = [];
            for b = 1:length(bonds)
                if (length(find(bonds(b).label)) <= delta)
                    badTracks = [badTracks,b];
                end
            end
 
            bonds(badTracks) = [];

            for b = 1:length(bonds)
                this.times{b} = find(bonds(b).label);
                this.label{b} = bonds(b).label(this.times{b});
                this.cells(b,:) = bonds(b).cells;
            end
            
        end
        
    end
               
end