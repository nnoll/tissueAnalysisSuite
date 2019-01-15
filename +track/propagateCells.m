function [ cTrack ] = propagateCells( cPair, T )
    % PROPAGATE CELLS 

    if (T(2) - T(1) == 1)
        cTrack = cPair{T(2)};
    else
        ind = 1:length(cPair{T(2)});
        cTrack = 1:length(cPair{T(2)});
        for tt = T(2):-1:(T(1)+1)
            cTrack(ind) = cPair{tt}(cTrack(ind));
            ind = find(cTrack);
        end
    end
    
end

