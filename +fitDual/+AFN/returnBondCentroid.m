function [ d0, rho ] = returnBondCentroid( Struct, i0 )
    % RETURN BONDS 

    d0 = [];
    rho = [];
    
    for b = 1:length(Struct.Bdat)
        if (all(ismember(Struct.Bdat(b).cells,i0)))
            
            row = zeros(1,length(i0));
            row(i0==Struct.Bdat(b).cells(1)) = 1;
            row(i0==Struct.Bdat(b).cells(2)) = -1;
            
            if (Struct.Bdat(b).radius < inf) % Finite curvature.                
                rBar = double(Struct.Bdat(b).rBar);
                
                if (size(rBar,1) == 1)
                    rBar = rBar';
                end
                
                rho = [rho,rBar];
                d0 = [d0;row];
                
            end
            
        end
    end   
    
    rho = rho';
    d0 = sparse(d0);
    
end

