function [ d0, t1, t2, r1, r2, blist, bLength ] = returnBonds( Struct, i0 )
    % RETURN BONDS 

    d0 = [];
    t1 = [];
    t2 = [];
    
    r1 = [];
    r2 = [];
    
    blist = [];
    bLength = [];
    for b = 1:length(Struct.Bdat)
        if (all(ismember(Struct.Bdat(b).cells,i0)))
            blist =  [blist,b];
            row = zeros(1,length(i0));
            row(i0==Struct.Bdat(b).cells(1)) = 1;
            row(i0==Struct.Bdat(b).cells(2)) = -1;
            d0 = [d0;row];
            
            rC = Struct.Cdat(Struct.Bdat(b).cells(1)).centroid.coord - Struct.Cdat(Struct.Bdat(b).cells(2)).centroid.coord;
            rB = double([Struct.Vdat(Struct.Bdat(b).verts(1)).vertxcoord;Struct.Vdat(Struct.Bdat(b).verts(1)).vertycoord]) - ...
                 double([Struct.Vdat(Struct.Bdat(b).verts(2)).vertxcoord;Struct.Vdat(Struct.Bdat(b).verts(2)).vertycoord]);
             
            if (sign(rC * [0,-1;1,0] * rB) < 0)
                v1 = Struct.Bdat(b).verts(2);
                v2 = Struct.Bdat(b).verts(1);
                rB = -rB;
            else
                v1 = Struct.Bdat(b).verts(1);
                v2 = Struct.Bdat(b).verts(2);
            end
            
            L = sqrt(sum(rB.^2));
            rB = rB/L;
            
            r1 = [r1,double([Struct.Vdat(v1).vertxcoord;Struct.Vdat(v1).vertycoord])];
            r2 = [r2,double([Struct.Vdat(v2).vertxcoord;Struct.Vdat(v2).vertycoord])];
            bLength = [bLength,L];

            if (Struct.Bdat(b).radius < inf) % Finite curvature.
                R = Struct.Bdat(b).radius;

                rBar = Struct.Bdat(b).rBar;
                if (size(rBar,1) == 1)
                    rBar = rBar';
                end
                
                tmp1 = r1(:,size(r1,2)) - rBar;
                tmp2 = r2(:,size(r2,2)) - rBar;
                if (sign(tmp1(1)*tmp2(2) - tmp1(2)*tmp2(1)) > 0)
                    t1 = [t1,[0,-1;1,0] * tmp1/(sqrt(sum(tmp1.^2)))];
                    t2 = [t2,[0,1;-1,0] * tmp2/(sqrt(sum(tmp2.^2)))];
                else
                    t1 = [t1,[0,1;-1,0] * tmp1/(sqrt(sum(tmp1.^2)))];
                    t2 = [t2,[0,-1;1,0] * tmp2/(sqrt(sum(tmp2.^2)))];
                end
                
            else % Flat edge
                t1 = [t1,-rB];
                t2 = [t2,rB];
            end
            
        end
    end   
    
    t1 = t1';
    t2 = t2';
    r1 = r1';
    r2 = r2';
    d0 = sparse(d0);
    
end

