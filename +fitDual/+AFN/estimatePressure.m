function [ p ] = estimatePressure( q, i0, Struct )
    % ESTIMATE PRESSURE 
    
    L = zeros(2*length(Struct.Bdat),size(q,1));
    badBonds = [];
    equBonds = [];
    nB = length(Struct.Bdat);

    for b = 1:nB
        if (all(ismember(Struct.Bdat(b).cells,i0)))
            ind1 = find(i0 == Struct.Bdat(b).cells(1));
            ind2 = find(i0 == Struct.Bdat(b).cells(2));
            if (Struct.Bdat(b).radius > 0)
                rho = Struct.Bdat(b).rBar;
%                 rhoN = sqrt(sum(rho.^2));
                L(b,ind1) = (q(ind1,1) - rho(1));
                L(b,ind2) = (rho(1) - q(ind2,1));
                L(b+nB,ind1) = (q(ind1,2) - rho(2));
                L(b+nB,ind2) = (rho(2) - q(ind2,2));
            else
                L(b,ind1) = 1;
                L(b,ind2) = -1;
                equBonds = [equBonds,b];
            end
        else
            badBonds = [badBonds,b];
        end
    end
    
    L([badBonds,nB+badBonds,nB+equBonds],:) = [];
    nE = size(L,1);
    L = [L;ones(1,size(L,2))];
   
    b = [zeros(nE,1);size(L,2)];
    p = L \ b;
        
end

