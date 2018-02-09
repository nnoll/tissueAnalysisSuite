function [ deltaR ] = fitMissingDisplacements( Struct, cTrack, T, Rq )
    % FIT MISSING DISPLACEMENTS 
    
    deltaRc = zeros(sum(cTrack>0),2);
    Rc = zeros(size(deltaRc));
    
    n = 1;
    badVerts = [];
    for c = 1:length(cTrack)
        if (cTrack(c) > 0 && c > 1 && cTrack(c) > 1)
            c1 = c;
            c2 = cTrack(c);
            deltaRc(n,:) = Struct(T(2)).Cdat(c1).centroid.coord(1:2) - Struct(T(1)).Cdat(c2).centroid.coord(1:2);
            Rc(n,:) = Struct(T(1)).Cdat(c2).centroid.coord(1:2);
            n = n + 1;
        elseif (cTrack(c)>0)
            badVerts = [badVerts,n];
            n = n + 1;
        end
    end

    deltaRc(badVerts,:) = [];
    Rc(badVerts,:) = [];
    
    Fx = scatteredInterpolant(Rc(:,1),Rc(:,2),deltaRc(:,1));
    Fy = scatteredInterpolant(Rc(:,1),Rc(:,2),deltaRc(:,2));
    
    deltaR = [Fx(Rq(:,1),Rq(:,2)),Fy(Rq(:,1),Rq(:,2))];
    
end

