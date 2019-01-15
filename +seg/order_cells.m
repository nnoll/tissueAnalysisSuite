function [ Struct ] = order_cells( Struct )
%ORDER_CELLS 

for t=1:length(Struct)
    for c=1:length(Struct(t).Cdat)
        if (~isempty(Struct(t).Cdat(c).cell_good) && Struct(t).Cdat(c).all_threefold && Struct(t).Cdat(c).cell_good == 1)
            ncells = Struct(t).Cdat(c).ncells;
            R = Struct(t).Cdat(c).centroid.coord;
            theta = zeros(1,length(ncells));
            for ii=1:length(ncells);
                Rc = Struct(t).Cdat(ncells(ii)).centroid.coord;
                delta = Rc - R;
                theta(ii) = atan2(delta(2),delta(1));
            end
            [~,ind] = sort(theta,'descend');
            Struct(t).Cdat(c).orderedncells = ncells(ind);
        end
    end

    for v=1:length(Struct(t).Vdat)
        if (length(Struct(t).Vdat(v).ncells)==3 && length(Struct(t).Vdat(v).nverts)==3 )
            ncells = Struct(t).Vdat(v).ncells;
            R = [Struct(t).Vdat(v).vertxcoord;Struct(t).Vdat(v).vertycoord];
            theta = zeros(1,3);
            for ii=1:length(ncells);
                Rc = Struct(t).Cdat(ncells(ii)).centroid.coord;
                delta = Rc' - R;
                theta(ii) = atan2(delta(2),delta(1));
            end
            [~,ind] = sort(theta,'descend');
            Struct(t).Vdat(v).orderedncells = ncells(ind);
        end
    end
end

end

