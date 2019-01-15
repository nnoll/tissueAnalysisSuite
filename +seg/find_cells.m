function [Cdat]=find_cells(Vdat,L)


%  S=regionprops(L,'Area','Orientation','Perimeter','Centroid','Eccentricity','PixelIdxList','MajorAxisLength','MinorAxisLength');

S=regionprops(L,'Centroid');

numberofcells=max(max(L));


% Cdat.ncells=[];
% Cdat.nverts=[];
% Cdat.nbonds=[];
% Cdat.area=[];
% Cdat.orientation=[];
% Cdat.perimeter=[];
% Cdat.centroid.coord=[];
% Cdat.eccentricity=[];
% Cdat.major=[];
% Cdat.minor=[];
% Cdat.Myopulsestats=[];
% Cdat.cellpixels=[];

Cdat(numberofcells) = struct('centroid',[],'nverts',[],'numv', 0);
for ii=1:numberofcells
    % Cdat(ii).area=S(ii).Area;
    % Cdat(ii).orientation=S(ii).Orientation;
    % Cdat(ii).perimeter=S(ii).Perimeter;
    Cdat(ii).centroid.coord = S(ii).Centroid;
    % Cdat(ii).eccentricity=S(ii).Eccentricity;
    % Cdat(ii).major=S(ii).MajorAxisLength;
    % Cdat(ii).minor=S(ii).MinorAxisLength;
    % Cdat(ii).cellpixels=S(ii).PixelIdxList;
    %%%% Cdat(ii).nverts=ceil(find([Vdat.ncells]==ii)./3);
    Cdat(ii).nverts = zeros(1,6);
    Cdat(ii).numv = 0;
end

for ii=1:length(Vdat)
    for jj=1:length(Vdat(ii).ncells)
        Cdat(Vdat(ii).ncells(jj)).nverts(Cdat(Vdat(ii).ncells(jj)).numv+1) = ii;
        Cdat(Vdat(ii).ncells(jj)).numv = Cdat(Vdat(ii).ncells(jj)).numv + 1;
    end
end

%%Delete zeros in vertex
for ii=1:numberofcells
    if (Cdat(ii).numv < 6 && Cdat(ii).numv > 0)
        Cdat(ii).nverts((Cdat(ii).numv+1):6) = [];
    elseif (Cdat(ii).numv == 0)
        Cdat(ii).nverts = [];
    end
end




    