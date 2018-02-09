function [ h ] = loadHessian( h, bCells, NC, Z )
%LOADHESSIAN Summary of this function goes here
%   Detailed explanation goes here

    pF = Z*(Z+1)/2;
    NB = size(bCells,1);
    row = zeros(pF*NB,1);
    col = zeros(pF*NB,1);
    val = zeros(pF*NB,1);
    n = 1;
    delta = Z/2 + 1;
    for ii = 1:Z
        for jj = ii:Z
            
            II = ((n-1)*NB) + (1:NB);
            row(II) = 1:NB;
            
            if (ii <= Z/2)
                row(II) = bCells(:,1) + (ii-1)*NC;
            else
                row(II) = bCells(:,2) + (ii-delta)*NC;
            end
            
            if (jj <= Z/2)
                col(II) = bCells(:,1) + (jj-1)*NC;
            else
                col(II) = bCells(:,2) + (jj-delta)*NC;
            end

            val(II) = h(:,ii,jj);
            
            n = n + 1;
        end
    end
    
    h = sparse(row(val~=0),col(val~=0),val(val~=0),Z/2*NC,Z/2*NC);
    
end

