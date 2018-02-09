function [ Ae ] = computeEdgeAreas( this )
%   COMPUTE EDGE AREAS 

    r = this.computePrimalVerts();
    Ae = zeros(size(this.d0,1),2);
    edge = this.returnEdges;
    for e = 1:size(edge,1)
        verts = find(this.d1(:,e)~=0);
        if (length(verts) == 2)
            delta1 = r(verts(1),:) - this.q(edge(e,1),:);
            delta2 = r(verts(2),:) - this.q(edge(e,1),:);
            
            Ae(e,1) = .5*abs(delta1(1)*delta2(2) - delta1(2)*delta2(1));
            
            delta1 = r(verts(1),:) - this.q(edge(e,2),:);
            delta2 = r(verts(2),:) - this.q(edge(e,2),:);
            
            Ae(e,2) = .5*abs(delta1(1)*delta2(2) - delta1(2)*delta2(1));
        end
    end
    
end

