function [ ] = inferCurve( Struct, mem )
%INFERCURVE Summary of this function goes here
%   Detailed explanation goes here

   for ii = 1:length(Struct.Bdat)
        clf
        ii
        imshow(mem,[])
        hold all
        viscircles(Struct.Bdat(ii).rBar',Struct.Bdat(ii).radius);
        scatter(Struct.Bdat(ii).rBar(1),Struct.Bdat(ii).rBar(2),'filled')
        plot([Struct.Vdat(Struct.Bdat(ii).verts).vertxcoord],[Struct.Vdat(Struct.Bdat(ii).verts).vertycoord],'LineWidth',2)
        pause
    end
    
end

