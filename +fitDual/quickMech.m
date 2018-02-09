function [ Struct, T, P, K, i0, i1 ] = quickMech( Struct, mode )
    % QUICK MECH 

    for t = 1:length(Struct)
        t
        if (mode == 1)
            [ ~, T{t}, ~, dV, iVerts, i0{t} ] = fitDual.ATN.fitTensionGraph( Struct(t), 1 );
            
            bVerts = zeros(length(Struct(t).Bdat),2);
            for b = 1:length(Struct(t).Bdat)
               bVerts(b,1:length(Struct(t).Bdat(b).verts)) = Struct(t).Bdat(b).verts;
            end
            
            i1{t} = zeros(length(T{t}),1);
            for b = 1:size(dV,1)
                verts = iVerts(dV(b,:)~=0);
                i1{t}(b) = find(ismember(bVerts,verts,'rows'));
            end
            
            P{t} = ones(length(i0{t}),1);
            K{t} = zeros(size(T{t}));
            
        else
            
            [ q, p, ~, bulk0, ext0 ] = fitDual.AFN.returnReducedDual( Struct(t), 1 );
            i0{t} = [bulk0,ext0];
            P{t} = p;
            [ T{t}, K{t}, i1{t} ] = pressure.estimateTensions( p, q, Struct(t) );
            
        end
    end
    
    for t = 1:length(Struct)
        
       for b = 1:length(i1{t})
          Struct(t).Bdat(i1{t}(b)).tension = T{t}(b);
       end
       
       for c = 1:length(i0{t})
           Struct(t).Cdat(i0{t}(c)).pressure = P{t}(c);
       end
       
    end
    
end

