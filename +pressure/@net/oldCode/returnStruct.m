function [ Struct ] = returnStruct( this, alpha )
    % RETURN STRUCT 
    
    if (nargin == 1 || alpha == 0)
        rV = this.computePrimalVerts();
        d0 = this.d0;
        d1 = this.d1;
        q = this.q;
        [ faces ] = this.computeFaces( rV );

        [ rho, R, edgeVerts, edgeCells ] = this.computeEdgeArcs(1);
    else
        rV = round(alpha*this.computePrimalVerts())/alpha;
        [ faces ] = this.computeFaces( this.computePrimalVerts() );
        [ rho, R, edgeVerts, edgeCells ] = this.computeEdgeArcs(1);
        
        d0 = this.d0;
        d1 = this.d1;
        
        contin = 1;
        D = pdist2(rV,rV);
        D = D + diag(ones(size(D,1),1)) + tril(ones(size(D,1)),-1);
        [row,col] = find(D==0);
        
        while (contin==1)
            delFaces = zeros(length(row),1);
            colFaces = zeros(length(row),1);
            delEdges = zeros(length(row),1);

            for ii = 1:length(row)
                edge = find( (d1(row(ii),:) .* d1(col(ii),:)) ~= 0 );
                if (~isempty(edge))
                    edgeFaces = find( d1(:,edge)~=0 );
                    colFaces(ii) = edgeFaces(1);
                    delFaces(ii) = edgeFaces(2);
                    delEdges(ii) = edge;
                end
            end

            ind = find(delFaces);
            [delFaces,sind] = unique(delFaces(ind));
            colFaces = colFaces(ind);
            colFaces = colFaces(sind);
            delEdges = delEdges(ind);
            delEdges = delEdges(sind);
            
            for ii = 1:length(colFaces)
                if (ismember(colFaces(ii),delFaces))
                   colFaces(ii) = colFaces((ismember(delFaces,colFaces(ii)))); 
                end
            end

            for ii = 1:length(colFaces)
                d1(colFaces(ii),:) = sum(d1([colFaces(ii),delFaces(ii)],:),1);
            end

            shiftFaces = zeros(size(faces));
            shiftEdgeVerts = zeros(size(edgeVerts));

            for ii = 1:length(colFaces)
               faces(faces == delFaces(ii)) = colFaces(ii); 
               edgeVerts(edgeVerts == delFaces(ii)) = colFaces(ii);
            end

            for ii = 1:length(colFaces)
                shiftFaces(faces > delFaces(ii)) = shiftFaces(faces > delFaces(ii)) - 1;
                shiftEdgeVerts(edgeVerts > delFaces(ii)) = shiftEdgeVerts(edgeVerts > delFaces(ii)) - 1; 
            end

            faces = faces + shiftFaces;
            edgeVerts = edgeVerts + shiftEdgeVerts;

            d1(delFaces,:) = [];
            d1(:,delEdges) = [];
%             D(delFaces,:) = [];
%             D(:,delFaces) = [];
            d0(delEdges,:) = [];
            rV(delFaces,:) = [];
            edgeVerts(delEdges,:) = [];
            edgeCells(delEdges,:) = [];
            R(delEdges) = [];
            rho(delEdges,:) = [];
            q = this.q;
            
            D = pdist2(rV,rV);
            D = D + diag(ones(size(D,1),1)) + tril(ones(size(D,1)),-1);
            [row,col] = find(D==0);
            
            if (isempty(row) || isempty(delFaces))
                contin = 0;
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
            Struct.Vdat(v).ncells = this.pTri(v,:);
        else
            conEdges = ( d1(v,:) ~= 0 );
            conCells = find( sum(abs(d0(conEdges,:)),1) ~= 0);
            deltaR = q(conCells,:);
            deltaR = bsxfun(@minus,deltaR,mean(deltaR,1));
            theta = atan2(deltaR(:,2),deltaR(:,1));
            theta = mod(theta,2*pi);
            [~,ind] = sort(theta);
            Struct.Vdat(v).ncells = conCells(ind);
        end
        Struct.Vdat(v).nverts = find(vertAdj(v,:));
    end
        
    for c = 1:size(cellAdj,1)
        Struct.Cdat(c).nverts = faces(c,~isnan(faces(c,:)));
        Struct.Cdat(c).centroid.coord = mean(rV(Struct.Cdat(c).nverts,:),1);
        Struct.Cdat(c).numverts = length(Struct.Cdat(c).nverts);
        Struct.Cdat(c).ncells = find(cellAdj(c,:));
    end
    
    for b = 1:size(d0,1)
        Struct.Bdat(b).cells = edgeCells(b,:);
        Struct.Bdat(b).verts = edgeVerts(b,:);

        if (R(b) < inf)
            Struct.Bdat(b).radius = R(b);
            Struct.Bdat(b).rBar = rho(b,:);
        else
            Struct.Bdat(b).radius = inf;
            Struct.Bdat(b).rBar = [inf,inf];
            
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
        
    end
    
    for c = 1:length(Struct.Cdat)
        Struct.Cdat(c).ncells = unique(Struct.Cdat(c).ncells);
        Struct.Cdat(c).nverts = unique(Struct.Cdat(c).nverts);
    end
        
end

