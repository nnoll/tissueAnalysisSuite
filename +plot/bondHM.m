function [ bondHM ] = bondHM( Struct, Lx, Ly, fieldName )
    % BOND HM 

    d = 1;
    maxInd = Ly + Ly*(Lx-1);
    
    bondHM = zeros(Ly,Lx,length(Struct));
    for t = 1:length(Struct)
        tmp = zeros(Ly,Lx);
        for b = 1:length(Struct(t).Bdat)
            v1 = Struct(t).Bdat(b).verts(1);
            v2 = Struct(t).Bdat(b).verts(2);
            if (v1 > 0 && v2 > 0)
                r1 = round(Struct(t).Vdat(v1).vertycoord) + Ly*(round(Struct(t).Vdat(v1).vertxcoord)-1);
                r2 = round(Struct(t).Vdat(v2).vertycoord) + Ly*(round(Struct(t).Vdat(v2).vertxcoord)-1);

                Pix = vertcat(Struct(t).Bdat(b).pix,[r1;r2]);
                nPix = length(Pix);
                dilatePix = zeros(1,(2*d+1)^2*nPix);
                k = 1;
                for x = -d:d
                    for y = -d:d
                        dilatePix(nPix*(k-1)+1:k*nPix) = Pix + y + Ly*x;
                        k = k + 1;
                    end
                end
 
                dilatePix = dilatePix( (dilatePix <= maxInd) .* (dilatePix > 0) == 1 );
                tmp(dilatePix) = min(Struct(t).Bdat(b).(fieldName),500);
            end
        end
        bondHM(:,:,t) = tmp;
    end


end

