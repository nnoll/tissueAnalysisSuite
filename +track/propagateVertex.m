function [ vTrack, cTrack ] = propagateVertex( vPair, cPair, Struct, T )
    % PROPAGATE VERTEX 
    
    if (T(2) - T(1) == 1)
        vTrack = vPair{T(2)};
        cTrack = cPair{T(2)};
    else
        % First propagate cells in time.
        cTrack = track.propagateCells(cPair,T);
        vTrack = zeros(length(Struct(T(2)).Vdat),1);
        
        % Use propagated cells to track vertices.
        zMax = 3;
        for v = 1:length(Struct(T(1)).Vdat)
            if (length(Struct(T(1)).Vdat(v).ncells) > zMax)
                zMax = length(Struct(T(1)).Vdat(v).ncells);
            end
        end
        
        for v = 1:length(Struct(T(2)).Vdat)
            if (length(Struct(T(2)).Vdat(v).ncells) > zMax)
                zMax = length(Struct(T(2)).Vdat(v).ncells);
            end
        end
        
        
        v0Cells = zeros(length(Struct(T(1)).Vdat),zMax);
        for v = 1:length(Struct(T(1)).Vdat)
            v0Cells(v,1:length(Struct(T(1)).Vdat(v).ncells)) = Struct(T(1)).Vdat(v).ncells;
        end
        
        vCells = zeros(length(Struct(T(2)).Vdat),zMax);
        for v = 1:length(Struct(T(2)).Vdat)
            vCells(v,1:length(Struct(T(2)).Vdat(v).ncells)) = cTrack(Struct(T(2)).Vdat(v).ncells);
        end
%         vCells = cTrack(vCells);
        
        v0Cells = sort(v0Cells,2);
        vCells = sort(vCells,2);
    
        [~,ind0,ind] = intersect(v0Cells,vCells,'rows');
        vTrack(ind) = ind0;
       
        badVerts1 = find(~ismember(1:length(Struct(T(1)).Vdat),vTrack));
        badVerts2 = find(vTrack == 0);
        
        R0 = [Struct(T(1)).Vdat.vertxcoord;Struct(T(1)).Vdat.vertycoord]';
        R = [Struct(T(2)).Vdat.vertxcoord;Struct(T(2)).Vdat.vertycoord]';

        Rb = R(badVerts2,:);
        Rb0 = R0(badVerts1,:);
        delta0 = isogonal.fitMissingDisplacements(Struct,cTrack,T,Rb0);
        
        D = pdist2(Rb,Rb0+delta0);
        D(D > 10) = inf;
        
        [mpair,~] = track.munkres(D);
        matchedVerts = find(mpair);
        
        vTrack(badVerts2(matchedVerts)) = badVerts1(mpair(matchedVerts));
        
    end

end

