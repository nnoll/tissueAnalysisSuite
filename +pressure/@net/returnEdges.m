function [ edge ] = returnEdges( this )
    % RETURN EDGES 

    edge = zeros(size(this.d0,1),2);
    for e = 1:size(this.d0,1)
        if (~isempty(this.cellLabels))
            edge(e,1) = this.cellLabels(find(this.d0(e,:)==1));
            edge(e,2) = this.cellLabels(find(this.d0(e,:)==-1));
        else
            edge(e,1) = find(this.d0(e,:)==1);
            edge(e,2) = find(this.d0(e,:)==-1);
        end
    end

    edge = sort(edge,2);
    
end

