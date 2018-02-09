function [ ] = inferCurveVertex( Struct, mem )
%INFERCURVE Summary of this function goes here
%   Detailed explanation goes here

   for v = 1:length(Struct.Vdat)
       if (length(Struct.Vdat(v).nverts) == 3 && all(Struct.Vdat(v).bond~=0)) % If it's three-fold
           clf
           imshow(mem,[])
           hold all
           v
           for b = Struct.Vdat(v).bond
                Struct.Bdat(b).rBar
                viscircles(Struct.Bdat(b).rBar',Struct.Bdat(b).radius);
                scatter(Struct.Bdat(b).rBar(1),Struct.Bdat(b).rBar(2),'r','filled')
                plot([Struct.Vdat(Struct.Bdat(b).verts).vertxcoord],[Struct.Vdat(Struct.Bdat(b).verts).vertycoord],'LineWidth',2,'Color','b')
           end
           scatter(Struct.Vdat(v).vertxcoord,Struct.Vdat(v).vertycoord,'g','filled')

           pause
       end
   end
    
end

