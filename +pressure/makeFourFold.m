function [ Struct ] = makeFourFold( Struct )
    % MAKE FOUR FOLD 

        rv(:,1) = [Struct.Vdat.vertxcoord];
        rv(:,2) = [Struct.Vdat.vertycoord];
        
        D = pdist2(rv,rv);
        
        D = D + diag(ones(size(D,1),1)) + tril(ones(size(D,1)),-1);
        [row,col] = find(D==0);
        
        if (~isempty(row))
            
            for ii = 1:length(row)
                Struct.Vdat(row(ii)).nverts = unique([Struct.Vdat(row(ii)).nverts,Struct.Vdat(col(ii)).nverts]);
                Struct.Vdat(row(ii)).ncells = unique([Struct.Vdat(row(ii)).ncells,Struct.Vdat(col(ii)).ncells]);
            end
            
            for c = 1:length(Struct.Cdat)
                if (any(ismember(col,Struct.Cdat(c).nverts)))
                    
                end
        end
end

