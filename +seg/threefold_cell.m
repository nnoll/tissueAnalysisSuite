function [ nStruct ] = threefold_cell( Struct )
%THREEFOLD_CELL Will append a boolean flag onto each cell that tells us if
%it contains only threefold vertices or not.

    nStruct = Struct;
    for t=1:length(Struct)
        for c = 1:length(Struct(t).Cdat)
            vert_list = Struct(t).Cdat(c).nverts;
            if (~isempty(vert_list)) 
                three_fold = length(vert_list) >= 1;
                for v=vert_list
                    three_fold = three_fold * (length(Struct(t).Vdat(v).nverts) == 3);
                end
                nStruct(t).Cdat(c).all_threefold = three_fold;
            else
                nStruct(t).Cdat(c).all_threefold = 0;
            end
        end
    end

end

