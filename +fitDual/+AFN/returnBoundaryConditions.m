function [ dBC, tBC, rBC ] = returnBoundaryConditions( Struct, bulk0, ext0, bulkVerts )
    % RETURN BOUNDARY CONDITIONS 

    i0 = [bulk0,ext0];
    dBC = [];
    tBC = [];
    rBC = [];
    
    for b = 1:length(Struct.Bdat)
        if (all(ismember(Struct.Bdat(b).cells,ext0))) %If bond is external.
            row = zeros(1,length(i0));
            row(i0==Struct.Bdat(b).cells(1)) = 1;
            row(i0==Struct.Bdat(b).cells(2)) = -1;
            dBC = [dBC;row];
            
            if (ismember(Struct.Bdat(b).verts(1),bulkVerts))
                v = Struct.Bdat(b).verts(1);
                vExt = Struct.Bdat(b).verts(2);
            else
                v = Struct.Bdat(b).verts(2);
                vExt = Struct.Bdat(b).verts(1);
            end
            
            rB = [Struct.Vdat(vExt).vertxcoord;Struct.Vdat(vExt).vertycoord] ... 
               - [Struct.Vdat(v).vertxcoord;Struct.Vdat(v).vertycoord];
            L = sqrt(sum(rB.^2));
            rB = rB/L;
            
            rBC = [rBC,double([Struct.Vdat(v).vertxcoord;Struct.Vdat(v).vertycoord])];

            if (Struct.Bdat(b).radius < inf) % Finite curvature.

                rBar = Struct.Bdat(b).rBar;
                
                if (size(rBar,1) == 1)
                    rBar = rBar';
                end
                
                tmp = rBC(:,size(rBC,2)) - rBar;
                tBC = [tBC,[0,-1;1,0] * tmp/(sqrt(sum(tmp.^2)))];
            else % Flat edge
                tBC = [tBC,rB];
            end
            
        end
    end   
    
    tBC = tBC';
    rBC = rBC';
    dBC = sparse(dBC);
    
end

