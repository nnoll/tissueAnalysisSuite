function [ Struct, q, theta ] = fakeStruct( N, alpha, beta, gamma )
    % FAKE STRUCT 

    if (nargin == 3)
        [ q, tri ] = generate.meshGrid( N, alpha );
        theta = abs(beta*randn(size(q,1),1));

        r = fitATN.returnVertexPositions(q,theta,tri);
    else
        [ q, tri ] = generate.meshGrid( N, alpha );
        theta = abs(beta*randn(size(q,1),1));
        p = abs(1 + gamma*randn(size(q,1),1));
        [pTri,tri] = fitAFN.returnPTri(q,p,tri);
        r = fitAFN.returnVertexPositions(q,theta,p,tri,pTri);
    end
    
    Mesh = MeshTriangle(tri,q,'Omega');
    [ d0, d1, b0 ] = generate.primal( Mesh );
    [ faces ] = generate.faces( r, d0, d1 );
    
    cellAdj = -d0'*d0; cellAdj = cellAdj - diag(diag(cellAdj));
    vertAdj = -d1*d1'; vertAdj = vertAdj - diag(diag(vertAdj));
    
    Struct = struct('Vdat',[],'Cdat',[]);
    
    % Initialize Vdat
    for v = 1:size(r,1)
       Struct.Vdat(v).vertxcoord = r(v,1);
       Struct.Vdat(v).vertycoord = r(v,2);
       Struct.Vdat(v).ncells = tri(v,:) + 1;
       Struct.Vdat(v).nverts = sort(find(vertAdj(v,:)));
    end
    
    % Initialize Cdat
    for c = 1:size(q,1)
        Struct.Cdat(c+1).ncells = sort(find(cellAdj(c,:)));
        Struct.Cdat(c+1).nverts = faces(c,~isnan(faces(c,:)));
        Struct.Cdat(c+1).centroid.coord = mean(r(Struct.Cdat(c+1).nverts,:),1)';
        
        Struct.Cdat(c+1).all_threefold = 1;
        for v = Struct.Cdat(c+1).nverts
            Struct.Cdat(c+1).all_threefold = Struct.Cdat(c+1).all_threefold * (length(Struct.Vdat(v).nverts) == 3);
        end
        
        if (ismember(c,b0))
           Struct.Cdat(c+1).ncells = [1,Struct.Cdat(c+1).ncells];
        end
    end
    Struct.Cdat(1).all_threefold = 0;
    Struct.Cdat(1).ncells = b0' + 1;
    Struct.Cdat(1).centroid.coord = [0,0];
end

