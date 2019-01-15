function [ bCells ] = bCells( d0, i0 )
    % B CELLS 

    bCells = zeros(size(d0,1),2);

    [bID,cPID] = find(d0==1);
    bCells(bID,1) = i0(cPID);

    [bID,cNID] = find(d0==-1);
    bCells(bID,2) = i0(cNID);

    bCells = sort(bCells,2);

end

