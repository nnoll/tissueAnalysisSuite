function [ Struct ] = bondChemistry( Struct, chem, d, name )
%RECORD_CHEMISTRY 

    maxInd = size(chem,1) + size(chem,1)*(size(chem,2)-1);

    if (~any(strcmp('chem',fieldnames(Struct(1).Bdat))))
        indX = 1;
    else
        indX = length(Struct(1).Bdat(1).name) + 1;
    end

    for t = 1:length(Struct)
        if (indX == 1)
            Struct(t).Bdat(1).name{1} = name;
        else
            Struct(t).Bdat(1).name{indX} = name;
        end
        chDat = squeeze(chem(:,:,t));
        for b = 1:length(Struct(t).Bdat)
            v1 = Struct(t).Bdat(b).verts(1);
            v2 = Struct(t).Bdat(b).verts(2);
            if (v1 > 0 && v2 > 0)
                r1 = round(Struct(t).Vdat(v1).vertycoord) + size(chem,1)*(round(Struct(t).Vdat(v1).vertxcoord)-1);
                r2 = round(Struct(t).Vdat(v2).vertycoord) + size(chem,1)*(round(Struct(t).Vdat(v2).vertxcoord)-1);

                Pix = vertcat(Struct(t).Bdat(b).pix,[r1;r2]);
                nPix = length(Pix);
                dilatePix = zeros(1,(2*d+1)^2*nPix);
                k = 1;
                for x = -d:d
                    for y = -d:d
                        dilatePix(nPix*(k-1)+1:k*nPix) = Pix + y + size(chem,1)*x;
                        k = k + 1;
                    end
                end
                
                dilatePix = dilatePix( (dilatePix <= maxInd) .* (dilatePix > 0) == 1 );
                Struct(t).Bdat(b).chem(indX) = mean(chDat(dilatePix));
            end
        end
    end

end

