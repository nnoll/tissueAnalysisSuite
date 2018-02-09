function [ Struct ] = bondLength( Struct )
%BONDLENGTH Summary of this function goes here
%   Detailed explanation goes here

    for t = 1:length(Struct)
        for b = 1:length(Struct(t).Bdat)
            v1 = Struct(t).Bdat(b).verts(1);
            v2 = Struct(t).Bdat(b).verts(2);
            Struct(t).Bdat(b).length = sqrt( sum( ...
            ([Struct(t).Vdat(v1).vertxcoord;Struct(t).Vdat(v1).vertycoord] - ...
             [Struct(t).Vdat(v2).vertxcoord;Struct(t).Vdat(v2).vertycoord]).^2) );
        end
    end

end

