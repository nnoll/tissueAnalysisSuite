function [ len, kappa ] = returnBondStatistics( Struct )
%RETURNBONDSTATISTICS 

    len = cell(length(Struct),1);
    kappa = cell(length(Struct),1);

    for t = 1:length(Struct)
        len{t} = zeros(length(Struct(t).Bdat),1);
        kappa{t} = zeros(length(Struct(t).Bdat),1);
        for b = 1:length(Struct(t).Bdat)
            v1 = Struct(t).Bdat(b).verts(1);
            v2 = Struct(t).Bdat(b).verts(2);
            
            r1 = [Struct(t).Vdat(v1).vertxcoord;Struct(t).Vdat(v1).vertycoord];
            r2 = [Struct(t).Vdat(v2).vertxcoord;Struct(t).Vdat(v2).vertycoord];
            
            len{t}(b) = sqrt(sum( (r1-r2).^2 ));
            
            kappa{t}(b) = len{t}(b)/Struct(t).Bdat(b).radius;
        end
    end

end

