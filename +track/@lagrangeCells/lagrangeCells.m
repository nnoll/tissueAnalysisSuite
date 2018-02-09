classdef lagrangeCells < handle
    %   Partitions bonds into unique identities to easily access their
    %   properties such as stored chemistry and lengths, etc.
    
    %---------------------------------------------------------------------
    % Properties
    %---------------------------------------------------------------------

    properties (SetAccess = protected)

        t0 % Initial time bond was observed
        label % Label in the data structure at the given time.
        delta % Cell must be tracked longer than delta to count.
        
    end


    %---------------------------------------------------------------------
    % Public methods
    %---------------------------------------------------------------------
    
    methods
        
        %------------------------------------------------------
        % Constructor
        %------------------------------------------------------
        
        function this = lagrangeCells( Struct, cPair, delta )
            [ cells ] = track.resolveUniqueCells( Struct, cPair );
            this.delta = delta;
            badTracks = [];
            for c = 1:length(cells)
                if (length(cells(c).label) <= delta)
                    badTracks = [badTracks,c];
                end
            end
            cells(badTracks) = [];
            
            this.t0 = [cells.t0];
            for c = 1:length(cells)
                this.label{c} = cells(c).label;
            end
        end
        
    end
               
end