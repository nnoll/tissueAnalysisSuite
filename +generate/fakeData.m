function [ Struct, L, r ] = fakeData( D, dim )
    % FAKE DATA 
    
    [ r ] = generate.triangularLattice( round(D), dim*[1,1] );
    r = r + 3*randn(size(r));
    d = abs(1*ones(size(r,2),1) + .15*randn(size(r,2),1));
    [ L ] = generate.weightedVoronoi( r, d, dim*[1,1] );
    [ L, Struct ] = seg.generate_structs(L,0,1);
    Struct = seg.threefold_cell(Struct);
    r = r';
    
end

