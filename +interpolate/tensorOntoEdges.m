function [ Struct ] = tensorOntoEdges( Struct, mTensor, sTensor )
%TENSORONTOEDGES 

    for t = 1:length(Struct)
        for b = 1:length(Struct(t).Bdat)
            if (~isempty(Struct(t).Bdat(b).tension))
               v1 = Struct(t).Bdat(b).verts(1);
               v2 = Struct(t).Bdat(b).verts(2);
               r1 = [Struct(t).Vdat(v1).vertxcoord,Struct(t).Vdat(v1).vertycoord];
               r2 = [Struct(t).Vdat(v2).vertxcoord,Struct(t).Vdat(v2).vertycoord];
               
               R = .5*(r1+r2);
               rB = r1 - r2;
               rB = rB / sqrt(sum(rB.^2));
               
               R = round(R/8);
               Struct(t).Bdat(b).myo = mTensor{t}(R(2),R(1),1) * rB(1) * rB(1) + ...
                                       mTensor{t}(R(2),R(1),2) * 2*rB(1) * rB(2) + ...
                                       mTensor{t}(R(2),R(1),3) * rB(2) * rB(2);
               Struct(t).Bdat(b).stress = sTensor{t}(R(2),R(1),1) * rB(1) * rB(1) + ...
                                       sTensor{t}(R(2),R(1),2) * 2*rB(1) * rB(2) + ...
                                       sTensor{t}(R(2),R(1),3) * rB(2) * rB(2);
            end
        end
    end


end

