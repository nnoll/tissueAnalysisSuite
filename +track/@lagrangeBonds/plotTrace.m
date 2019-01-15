function [ trace ] = plotTrace( this, id, cells, L, Struct )
%PLOTTRACE 

    trace = permute(repmat(L==0,[1,1,1,3]),[1,2,4,3]);
    
    ID1 = this.cells(id,1);
    ID2 = this.cells(id,2);
    
    d = 2;
    
    for t = 1:size(L,3)
        red = zeros(size(L,1),size(L,2));
        green = zeros(size(L,1),size(L,2));
        blue = zeros(size(L,1),size(L,2));
        
        if (t >= cells.t0(ID1) && t <= (cells.t0(ID1) + length(cells.label{ID1}) - 1))
            jj = t - cells.t0(ID1) + 1;
            red(L(:,:,t) == cells.label{ID1}(jj)) = 1;
        end
        
        if (t >= cells.t0(ID2) && t <= (cells.t0(ID2) + length(cells.label{ID2}) - 1))
            jj = t - cells.t0(ID2) + 1;
            green(L(:,:,t) == cells.label{ID2}(jj)) = 1;
        end
        
        tFind = (this.times{id} == t);
        if (sum(tFind) == 1)
            Pix = Struct(t).Bdat(this.label{id}(tFind)).pix;
            nPix = length(Pix);
            dilatePix = zeros(1,(2*d+1)^2*nPix);
            k = 1;
            for x = -d:d
                for y = -d:d
                    dilatePix(nPix*(k-1)+1:k*nPix) = Pix + y + size(blue,1)*x;
                    k = k + 1;
                end
            end
            blue(dilatePix) = 1;
        end
        
        trace(:,:,1,t) = (1-blue).*(trace(:,:,1,t) + red);
        trace(:,:,2,t) = (1-blue).*(trace(:,:,2,t) + green);
        trace(:,:,3,t) = trace(:,:,3,t) + blue;
        
    end
end

