function [ Struct ] = returnStruct( this, alpha )
    % RETURN STRUCT 
    
    if (nargin == 1 || alpha == 0)
        rV = this.computePrimalVerts();
        d0 = this.d0;
        d1 = this.d1;
        q = this.q;
        [ faces ] = this.computeFaces( rV );
        T = this.returnTension();
        [ rho, R, edgeVerts, edgeCells ] = this.computeEdgeArcs(1);
    else
        rV = this.computePrimalVerts();
        d0 = this.d0;
        d1 = this.d1;
        q = this.q;
        T = this.returnTension();
        [ faces ] = this.computeFaces( rV );
        rV = rV + alpha*randn(size(rV));

        [ rho, R, edgeVerts, edgeCells ] = this.computeEdgeArcs(1);
        rho = rho + alpha*randn(size(rho));
        for b = 1:length(R)
           if (R(b) < inf)
               R(b) = .5*( sqrt( sum( (rV(edgeVerts(b,1),:) - rho(b,:)).^2 ,2) ) + sqrt( sum( (rV(edgeVerts(b,1),:) - rho(b,:)).^2 ,2) ) );
           end
        end
    end
    
    cellAdj = d0'*d0;
    cellAdj = cellAdj - diag(diag(cellAdj));
    cellAdj(cellAdj ~= 0) = 1;

    vertAdj = d1*d1';
    vertAdj = vertAdj - diag(diag(vertAdj));
    vertAdj(vertAdj ~= 0) = 1; 
    
    Struct = struct('Vdat',[],'Cdat',[],'Bdat',[]);
    
    for v = 1:size(rV,1)
        Struct.Vdat(v).vertxcoord = rV(v,1);
        Struct.Vdat(v).vertycoord = rV(v,2);
        if (nargin == 1 || alpha == 0)
            Struct.Vdat(v).ncells = this.pTri(v,:) + 1;
        else
            conEdges = ( d1(v,:) ~= 0 );
            conCells = find( sum(abs(d0(conEdges,:)),1) ~= 0);
            deltaR = q(conCells,:);
            deltaR = bsxfun(@minus,deltaR,mean(deltaR,1));
            theta = atan2(deltaR(:,2),deltaR(:,1));
            theta = mod(theta,2*pi);
            [~,ind] = sort(theta);
            Struct.Vdat(v).ncells = conCells(ind) + 1;
        end
        Struct.Vdat(v).nverts = find(vertAdj(v,:));
    end
        
    Struct.Cdat(1).centroid.coord = [0,0];

    for c = 1:size(cellAdj,1)
        Struct.Cdat(c+1).nverts = faces(c,~isnan(faces(c,:)));
        Struct.Cdat(c+1).centroid.coord = mean(rV(Struct.Cdat(c+1).nverts,:),1);
        Struct.Cdat(c+1).numverts = length(Struct.Cdat(c+1).nverts);
        Struct.Cdat(c+1).ncells = find(cellAdj(c,:)) + 1;
        Struct.Cdat(c+1).actual_pressure = this.p(c);
        Struct.Cdat(c+1).q = q(c,:);
        Struct.Cdat(c+1).theta = this.theta(c);
    end
    
    for b = 1:size(d0,1)
        Struct.Bdat(b).cells = edgeCells(b,:) + 1;
        Struct.Bdat(b).verts = edgeVerts(b,:);
        Struct.Bdat(b).actual_tension = T(b);
        if (R(b) < inf)
            Struct.Bdat(b).radius = R(b);
            Struct.Bdat(b).rBar = rho(b,:);
        else
            Struct.Bdat(b).radius = inf;
            Struct.Bdat(b).rBar = [inf,inf];
        end
        
        % Need to create a new vertex.
        if (length(find(Struct.Bdat(b).verts)) < 2)
            NV = length(Struct.Vdat);
            Struct.Bdat(b).verts(2) = NV + 1;
            deltaQ = diff(q(edgeCells(b,:),:),1);
            deltaQ = .5*([0,-1;1,0]*deltaQ')/sqrt(sum(deltaQ.^2));

            Struct.Vdat(NV+1).vertxcoord = rV(Struct.Bdat(b).verts(1),1) + deltaQ(1);
            Struct.Vdat(NV+1).vertycoord = rV(Struct.Bdat(b).verts(1),2) + deltaQ(2);
            if ((Struct.Vdat(NV+1).vertxcoord^2 + Struct.Vdat(NV+1).vertycoord^2) <= sum(rV(Struct.Bdat(b).verts(1),:).^2,2))
                Struct.Vdat(NV+1).vertxcoord = rV(Struct.Bdat(b).verts(1),1) - deltaQ(1);
                Struct.Vdat(NV+1).vertycoord = rV(Struct.Bdat(b).verts(1),2) - deltaQ(2);
            end

            Struct.Vdat(NV+1).ncells = Struct.Bdat(b).cells;
            Struct.Vdat(NV+1).nverts = Struct.Bdat(b).verts(1);
            Struct.Vdat(Struct.Bdat(b).verts(1)).nverts = [Struct.Vdat(Struct.Bdat(b).verts(1)).nverts,NV+1];
            Struct.Cdat(Struct.Bdat(b).cells(1)).nverts = unique([Struct.Cdat(Struct.Bdat(b).cells(1)).nverts,NV+1]);
            Struct.Cdat(Struct.Bdat(b).cells(2)).nverts = unique([Struct.Cdat(Struct.Bdat(b).cells(2)).nverts,NV+1]);
        end
    end
    
    for c = 1:length(Struct.Cdat)
        Struct.Cdat(c).ncells = unique(Struct.Cdat(c).ncells);
        Struct.Cdat(c).nverts = unique(Struct.Cdat(c).nverts);
        Struct.Cdat(c).centroid.coord = [mean([Struct.Vdat(Struct.Cdat(c).nverts).vertxcoord]),mean([Struct.Vdat(Struct.Cdat(c).nverts).vertycoord])];
    end
        
end

