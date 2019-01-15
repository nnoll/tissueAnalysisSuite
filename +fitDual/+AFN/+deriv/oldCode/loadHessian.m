function [ h ] = loadHessian( h, bCells, NC )
%LOADHESSIAN Summary of this function goes here
%   Detailed explanation goes here

    NB = size(bCells,1);
    row = zeros(36*NB,1);
    col = zeros(36*NB,1);
    val = zeros(36*NB,1);
    n = 1;
    for ii = 1:8
        for jj = ii:8
            
            II = ((n-1)*NB) + (1:NB);
            row(II) = 1:NB;
            
            if (ii <= 4)
                rowH = bCells(:,1) + (ii-1)*NC;
            else
                rowH = bCells(:,2) + (ii-5)*NC;
            end
            
            if (jj <= 4)
                colH = bCells(:,1) + (jj-1)*NC;
            else
                colH = bCells(:,2) + (jj-5)*NC;
            end
            
            col(II) = rowH + (4*NC)*(colH-1);
            val(II) = h(:,ii,jj);
            
            n = n + 1;
        end
    end
    h = sparse(row(val~=0),col(val~=0),val(val~=0),NB,(4*NC)^2);
    
end

