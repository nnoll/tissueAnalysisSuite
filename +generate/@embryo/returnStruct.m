function [ Struct ] = returnStruct( this, alpha )
    % RETURN STRUCT 
    
    if (nargin == 1 || alpha == 0)
        rV = this.computeCellVerts();
        d0 = this.d0;
        d1 = this.d1;
        q = this.Mesh.Points;
        T = d0 * q; 
        T = sqrt( sum(T.^2,2) );

        badBonds = abs(d1'*rV(:,1)) > pi/2; % Introduce a cut in the tangent plane to make graph planar.
        badCells = sum(abs(d0(badBonds,:)),1) >= 1;
        
        T(badBonds) = [];
        d0(badBonds,:) = [];
        d1(:,badBonds) = [];
        
        d0(:,badCells) = [];
        q(badCells,:) = [];

        badVerts = (sum(abs(d1),2) == 1);
        singleBonds = sum(d1(badVerts,:),1) >= 1;
        rV(badVerts,:) = [];
        d1(badVerts,:) = [];
        d1(:,singleBonds) = [];
        d0(singleBonds,:) = [];
        T(singleBonds) = [];
        
        badBonds = sum(abs(d0),2) == 0;
        d0(badBonds,:) = [];
        d1(:,badBonds) = [];
        T(badBonds) = [];

        [ edgeCells, edgeVerts ] = this.returnEdges( d0, d1 );
        [ faces ] = this.computeFaces( rV, d0, d1 );
    else
        rV = this.computePrimalVerts();
        d0 = this.d0;
        d1 = this.d1;
        T = d0 * this.Mesh.Points; 
        T = sqrt( sum(T.^2,2) );
        [ faces, badFaces ] = this.computeFaces( rV );
        rV = rV + alpha*randn(size(rV));
    end
    
    cellAdj = d0'*d0;
    cellAdj = cellAdj - diag(diag(cellAdj));
    cellAdj(cellAdj ~= 0) = 1;

    vertAdj = d1*d1';
    vertAdj = vertAdj - diag(diag(vertAdj));
    vertAdj(vertAdj ~= 0) = 1; 
    
    Struct = struct('Vdat',[],'Cdat',[],'Bdat',[]);
    
    Struct.Cdat(1).nverts = [];
    Struct.Cdat(1).ncells = [];

    for v = 1:size(rV,1)
        Struct.Vdat(v).vertxcoord = rV(v,1);
        Struct.Vdat(v).vertycoord = rV(v,2);

        conEdges = ( d1(v,:) ~= 0 );
        conCells = find( sum(abs(d0(conEdges,:)),1) ~= 0);
        deltaR = q(conCells,:);
        deltaR = bsxfun(@minus,deltaR,mean(deltaR,1));
        theta = atan2(deltaR(:,2),deltaR(:,1));
        theta = mod(theta,2*pi);
        [~,ind] = sort(theta);
        Struct.Vdat(v).ncells = conCells(ind) + 1;
        if (length(Struct.Vdat(v).ncells) < 3)
            Struct.Vdat(v).ncells = [1,Struct.Vdat(v).ncells];
            Struct.Cdat(1).nverts = [Struct.Cdat(1).nverts,v];
            Struct.Cdat(1).ncells = [Struct.Cdat(1).ncells,Struct.Vdat(v).ncells(Struct.Vdat(v).ncells~=1)];
%             this.plotDual()
%             hold on
%             scatter(rV(v,1),rV(v,2),'r')
%             pause
        end
        Struct.Vdat(v).nverts = find(vertAdj(v,:));
    end
        
    Struct.Cdat(1).centroid.coord = [0,0];

    for c = 1:size(cellAdj,1)
        Struct.Cdat(c+1).nverts = faces(c,~isnan(faces(c,:)));
        Struct.Cdat(c+1).centroid.coord = mean(rV(Struct.Cdat(c+1).nverts,:),1);
        Struct.Cdat(c+1).numverts = length(Struct.Cdat(c+1).nverts);
        Struct.Cdat(c+1).ncells = find(cellAdj(c,:)) + 1;
        if (ismember(c+1,Struct.Cdat(1).ncells) && ~ismember(1,Struct.Cdat(c+1).ncells))
            Struct.Cdat(c+1).ncells = [1,Struct.Cdat(c+1).ncells];
        end
        Struct.Cdat(c+1).pressure = 1;
    end
    
    for b = 1:size(d0,1)
        Struct.Bdat(b).cells = edgeCells(b,:) + 1;
        if (ismember(0,edgeCells(b,:)))
           bndryCell = edgeCells(b,edgeCells(b,:)>0);
           Struct.Cdat(1).ncells = [Struct.Cdat(1).ncells,bndryCell+1];
        end
        Struct.Bdat(b).verts = edgeVerts(b,:);
        Struct.Bdat(b).actual_tension = T(b);

        Struct.Bdat(b).radius = inf;
        Struct.Bdat(b).rBar = [inf,inf];  
    end
    
    for c = 1:length(Struct.Cdat)
        Struct.Cdat(c).ncells = unique(Struct.Cdat(c).ncells);
        Struct.Cdat(c).nverts = unique(Struct.Cdat(c).nverts);
    end
    
    Struct = seg.threefold_cell(Struct);
    Struct.embed = @(az,el)[sin(el).*cos(az),sin(el).*sin(az),cos(el)];
end

