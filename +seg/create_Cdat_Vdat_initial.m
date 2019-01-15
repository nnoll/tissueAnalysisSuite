function [Struct]=create_Cdat_Vdat_initial(Ltrack,mode,threefold, very_far)

    % 
    %
    % Returns
    % -------
    
    for ii=1:length(Ltrack(1,1,:))
        clear L
        L=Ltrack(:,:,ii);
        [Vdat]= seg.find_vertices(L,threefold, very_far);
        [Cdat]=seg.find_cells(Vdat,L);
        if (mode==1)
            [Vdat]=seg.find_bonds(Vdat,L);
        end 
        Struct(ii).Vdat = Vdat;
        Struct(ii).Cdat = Cdat;
        clear Vdat Cdat
    end

    for ii=1:length(Struct)
        for jj=1:length(Struct(ii).Cdat) 
            Struct(ii).Cdat(jj).ncells=setdiff(unique([Struct(ii).Vdat(Struct(ii).Cdat(jj).nverts).ncells]),jj);
        end
    end

    % [Struct(ii)]=seg.findCneigh(Struct(ii));

end




