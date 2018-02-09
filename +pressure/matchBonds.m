function [ bInd ] = matchBonds( PN, bCells )
    % MATCH BONDS 
    
    edge = PN.returnEdges();

    [~,~,tmp] = intersect(edge,bCells,'rows'); 
    bInd(ismember(edge,bCells,'rows')) = tmp;
end

