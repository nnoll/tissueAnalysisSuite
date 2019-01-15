function [ vertexPairs, cellPairs ] = vertex( Struct, L, cutoff )
%Will `align' vertices over succesive/different time points - i.e. will 
%propagate labels of vertices in time.

    %Inputs:
    %1. Struct = Data structure containing observables at all times 
    %2. L = Label matrix from watershed.
    %2. cutoff = Length scale over which we won't ascribe a vertex identity 

    %Outputs:
    %1. vertexPairs = Tracked vertex identities. 

    %Compute cell pairs to be used in vertex mapping problem.
    [ cellPairs ] = track.cells( L, cutoff );
    
    %Initialize vertex pair data structure. 
    vertexPairs = cell(size(cellPairs));
    for t = 2:length(vertexPairs)
        % First pair vertices if you can find all tracked cells with the
        % same vertex in neighboring time-points.
        vertexPairs{t} = zeros(length(Struct(t).Vdat),1);
        
        % Repeat with current vertices but mapped to previous time-point.
        for v = 1:length(Struct(t).Vdat)
           nCells = cellPairs{t}(Struct(t).Vdat(v).ncells)'; 
           for vv = 1:length(Struct(t-1).Vdat)
                if ((length(nCells) == length(Struct(t-1).Vdat(vv).ncells)) && all(nCells == Struct(t-1).Vdat(vv).ncells))
                    vertexPairs{t}(v) = vv;
                end
           end
        end
        
        % Match any unpaired current vertices with unpaired old vertices.
        unpairedCurrentVerts = find(vertexPairs{t}==0);
        unpairedOldVerts = 1:length(Struct(t-1).Vdat);
        unpairedOldVerts = unpairedOldVerts(~ismember(unpairedOldVerts,vertexPairs{t}));
        
        clear R R0
        R(:,1) = [Struct(t).Vdat.vertxcoord]; R(:,2) = [Struct(t).Vdat.vertycoord];
        R = R(unpairedCurrentVerts,:);
        
        R0(:,1) = [Struct(t-1).Vdat.vertxcoord]; R0(:,2) = [Struct(t-1).Vdat.vertycoord];
        R0 = R0(unpairedOldVerts,:);
        
        D = pdist2(double(R),double(R0));
        D(D > cutoff) = inf;
        [mpair,~] = track.munkres(D);
        
        matchedVerts = find(mpair);
        vertexPairs{t}(unpairedCurrentVerts(matchedVerts)) = unpairedOldVerts(mpair(matchedVerts))';
    end


end